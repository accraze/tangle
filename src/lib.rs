//! A spectral smear for acoustic instruments.
//!
//! Signal chain:
//! 1. STFT: `FFT_SIZE` frames every `HOP` samples, periodic sqrt-Hann windows on analysis and
//!    synthesis, so overlap-add reconstructs the input exactly. nih-plug's `StftHelper` does the
//!    buffering; its latency is one frame.
//! 2. Once per frame, the knobs are read and the circuit's unitary U is rebuilt (`tangle_core`)
//!    from the seeded base angles, offset by Twist and Drift, entanglers scaled by Coupling, all
//!    scaled by Amount.
//! 3. Bins from 1 up to `TOP_BIN` are split into blocks of 2^qubits contiguous bands of K bins
//!    (`band_bins`, set by Width; the few bins left over at the top stay put). In each block, for
//!    each offset k, the vector of bin k from every band is multiplied by U. U is unitary, so this
//!    keeps the frame's energy. Mix blends towards the identity in the spectrum, which equals a dry
//!    blend; Transient Preserve pulls that blend down on each channel's attacks (`Detector`).
//!    Octaves mode replaces this step with `shift_octaves`: U's index is the octave instead.
//! 4. Inverse FFT and overlap-add.

use nih_plug::prelude::*;
use nih_plug_egui::widgets::generic_ui::{self, GenericSlider};
use nih_plug_egui::{create_egui_editor, egui, EguiState};
use realfft::num_complex::Complex32;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::f32::consts::{PI, TAU};
use std::sync::Arc;
use tangle_core::{Angles, Unitary, MAX_DIM, MAX_LAYERS, MAX_QUBITS};

const FFT_SIZE: usize = 2048;
const OVERLAP: usize = 4;
const HOP: usize = FFT_SIZE / OVERLAP;
/// Highest bin that can be mixed. DC (0) and Nyquist (FFT_SIZE / 2) must stay real for the
/// inverse FFT, so they're never mixed.
const TOP_BIN: usize = FFT_SIZE / 2 - 1;
/// Seed for the circuit's base angles. Changing it changes the sound of saved projects.
const SEED: u64 = 1;

/// Transient detection: rise in high-frequency content (dB per hop) where ducking starts, the
/// further rise that ducks fully, and how long the duck takes to let go (seconds, 1/e).
const ONSET_DB: f32 = 3.0;
const ONSET_RANGE_DB: f32 = 6.0;
const TRANSIENT_RELEASE: f32 = 0.1;
/// High-frequency content below this (roughly -80 dBFS) counts as silence, so the noise floor
/// doesn't trigger
const HFC_FLOOR: f32 = 1.0;

/// Octaves mode: 3 qubits index 8 octaves, bins OCTAVE_BASE·2^b to OCTAVE_BASE·2^(b+1), which
/// end exactly at TOP_BIN. That's every octave a 2048-point FFT can resolve, so the Qubits knob
/// doesn't apply. Bins below OCTAVE_BASE (under ~94 Hz at 48 kHz) stay put.
const OCTAVE_QUBITS: usize = 3;
const OCTAVE_BASE: usize = 4;

#[derive(Enum, PartialEq, Clone, Copy)]
enum Mode {
    Bands,
    Octaves,
}

/// Knob values for one STFT frame
#[derive(Clone, Copy, Default)]
struct Settings {
    amount: f32,
    mix: f32,
    /// 0 = 4-bin bands repeated up the spectrum, 1 = bands as wide as fit
    width: f32,
    qubits: usize,
    layers: usize,
    /// -1..1, offsets every rotation by up to ±π
    twist: f32,
    /// Entangler strength, 1 = full CNOTs
    coupling: f32,
    /// Drift depth: each rotation wobbles by up to ±π
    drift: f32,
    /// Hz
    drift_rate: f32,
    transient: f32,
    /// 0 = both channels share one circuit, 1 = the right channel runs an unrelated one
    spread: f32,
    /// Octaves mode instead of Bands
    octaves: bool,
    /// Octaves mode: furthest a partial can move, in octaves (1..=7)
    span: usize,
}

impl Settings {
    /// The circuit's angles at `drift_phase` (in cycles of the drift rate). Everything is scaled
    /// by Amount, so Amount 0 is exactly the identity whatever the other knobs say.
    fn angles(&self, base: &Angles, drift_phase: f64) -> Angles {
        let mut a = *base;
        for l in 0..MAX_LAYERS {
            for q in 0..MAX_QUBITS {
                let i = l * MAX_QUBITS + q;
                let offset = |i| (self.twist + self.drift * lfo(i, drift_phase)) * PI;
                a.ry[l][q] = (a.ry[l][q] + offset(i)) * self.amount;
                a.rz[l][q] = (a.rz[l][q] + offset(i + MAX_LAYERS * MAX_QUBITS)) * self.amount;
                a.ent[l][q] *= self.coupling * self.amount;
            }
        }
        a
    }

    /// Bins per band. Always a multiple of OVERLAP, so moving a bin by whole bands shifts it by a
    /// multiple of FFT_SIZE / HOP bins: its phase still advances consistently from hop to hop,
    /// overlapping frames add up coherently, and the move is a clean frequency shift. Otherwise
    /// the frames partly cancel (about -4 dB at full Amount) and the result gets phasey.
    fn band_bins(&self, dim: usize) -> usize {
        let max = TOP_BIN / dim / OVERLAP;
        let units = (max as f32).powf(self.width).round() as usize;
        units.clamp(1, max) * OVERLAP
    }
}

/// `a` moved toward `b` by `t`, so Spread morphs smoothly from one circuit to the other
fn blend(a: &Angles, b: &Angles, t: f32) -> Angles {
    let mut out = *a;
    for l in 0..MAX_LAYERS {
        for q in 0..MAX_QUBITS {
            out.ry[l][q] += (b.ry[l][q] - a.ry[l][q]) * t;
            out.rz[l][q] += (b.rz[l][q] - a.rz[l][q]) * t;
            out.ent[l][q] += (b.ent[l][q] - a.ent[l][q]) * t;
        }
    }
    out
}

/// A sine LFO per gate, each at its own rate (0.5x to 1.5x the drift rate) and phase, spread by
/// irrational steps so they never line up.
fn lfo(gate: usize, phase: f64) -> f32 {
    let rate = 0.5 + (gate as f64 * 0.618_033_988_75).fract();
    let offset = (gate as f64 * std::f64::consts::SQRT_2).fract();
    (std::f64::consts::TAU * (phase * rate + offset)).sin() as f32
}

/// Multiplies each cross-band vector of `bins` by `u`, blended with the input by `mix`. `bins` is
/// a whole number of blocks of `u.dim` bands of `k_len` bins.
fn mix_bands(bins: &mut [Complex32], u: &Unitary, k_len: usize, mix: f32) {
    for block in bins.chunks_exact_mut(u.dim * k_len) {
        for k in 0..k_len {
            let mut x = [Complex32::new(0.0, 0.0); MAX_DIM];
            for (b, x) in x.iter_mut().enumerate().take(u.dim) {
                *x = block[b * k_len + k];
            }
            for r in 0..u.dim {
                let ux: Complex32 = (0..u.dim).map(|c| u.cols[c][r] * x[c]).sum();
                block[r * k_len + k] = x[r] + (ux - x[r]) * mix;
            }
        }
    }
}

/// Per-channel phase vocoder state for Octaves mode, one entry per bin
#[derive(Default)]
struct Vocoder {
    /// Each bin's phase in the previous frame, to measure its frequency
    prev_phase: Vec<f32>,
    /// Output phase of the moved partial whose peak landed on each bin: last frame's, and this
    /// frame's being written
    out_phase: Vec<f32>,
    next_phase: Vec<f32>,
    /// This frame's moved partials, and how much of each bin stays put
    moved: Vec<Complex32>,
    /// Amplitude of the loudest partial landing on each bin this frame. It alone advances that
    /// bin's `out_phase`: a quieter partial (or a noise-floor peak) rounding to the same bin would
    /// otherwise advance it at its own rate and detune the loud one.
    loudest: Vec<f32>,
    stay: Vec<f32>,
    peaks: Vec<usize>,
    /// `window_response` sampled every 1/WINDOW_STEPS bin over ±(LOBE + 1) bins
    window: Vec<Complex32>,
}

impl Vocoder {
    fn new() -> Self {
        let bins = FFT_SIZE / 2 + 1;
        Self {
            prev_phase: vec![0.0; bins],
            out_phase: vec![0.0; bins],
            next_phase: vec![0.0; bins],
            moved: vec![Complex32::new(0.0, 0.0); bins],
            loudest: vec![0.0; bins],
            stay: vec![1.0; bins],
            peaks: Vec::with_capacity(bins),
            window: (0..=2 * (LOBE + 1) * WINDOW_STEPS)
                .map(|i| window_response(i as f32 / WINDOW_STEPS as f32 - (LOBE + 1) as f32))
                .collect(),
        }
    }

    fn reset(&mut self) {
        self.prev_phase.fill(0.0);
        self.out_phase.fill(0.0);
    }
}

/// `window_response(nu)` from `Vocoder::window`, for |nu| < LOBE + 1. The closed form is much
/// slower, and this runs 9 times per moved partial per octave.
fn window_at(table: &[Complex32], nu: f32) -> Complex32 {
    let x = (nu + (LOBE + 1) as f32) * WINDOW_STEPS as f32;
    let i = x as usize;
    table[i] + (table[i + 1] - table[i]) * (x - i as f32)
}

/// Bins either side of a redrawn partial's centre. The sine window's main lobe is ±1.5 bins and
/// its sidelobes fall off fast, so ±4 holds over 99.9% of the energy.
const LOBE: usize = 4;
/// Table resolution for `Vocoder::window_at`. Linear interpolation errs by about 1e-5.
const WINDOW_STEPS: usize = 256;
/// Energy of a redrawn unit partial's bins: N·Σw² by Parseval, a hair over what ±LOBE bins hold
const WINDOW_ENERGY: f32 = (FFT_SIZE * FFT_SIZE / 2) as f32;

/// The sqrt-Hann (sine) window's spectrum `nu` bins away from a partial, for a partial of
/// amplitude 1 and phase 0 at the frame's first sample
fn window_response(nu: f32) -> Complex32 {
    use realfft::num_complex::Complex64;
    let n = FFT_SIZE as f64;
    // Sum of e^(i 2π a j / N) over j in 0..N
    let sum = |a: f64| {
        let den = Complex64::from_polar(1.0, std::f64::consts::TAU * a / n) - 1.0;
        if den.norm() < 1e-12 {
            Complex64::new(n, 0.0)
        } else {
            (Complex64::from_polar(1.0, std::f64::consts::TAU * a) - 1.0) / den
        }
    };
    let nu = nu as f64;
    // sin(πj/N) = (e^(iπj/N) - e^(-iπj/N)) / 2i
    let w = (sum(0.5 - nu) - sum(-0.5 - nu)) / Complex64::new(0.0, 2.0);
    Complex32::new(w.re as f32, w.im as f32)
}

fn octave(bin: usize) -> usize {
    (bin / OCTAVE_BASE).ilog2() as usize
}

/// Octaves mode. U's index is the octave, and |U[r][c]|² is the share of each partial's energy
/// that moves from octave c to octave r (the Born rule: U is unitary, so each column's shares sum
/// to 1). The share that stays keeps its original bins, so U = I is exactly the identity.
///
/// Each spectral peak is a partial: the phase vocoder measures its true frequency, and dividing the
/// peak bin by the window's response there gives its amplitude. Each moved copy is redrawn as the
/// window's response centred on exactly 2^(r-c) times that frequency, with its phase advancing at
/// that rate from frame to frame, so copies are in tune and keep their level. Whole-bin shifts
/// can't do that: the copy's in-frame frequency stays off by the rounding, which detunes it.
///
/// shortcut: only peaks move, so the moved share of noisy, non-peak content (pick noise, breath)
/// is dropped rather than moved; split off a residual and shift it too if that's audible.
///
/// Only octaves within `span` of a partial's own get a share. The kept shares are scaled up so they
/// still sum to 1, which keeps the energy and leaves the circuit's proportions between them intact.
fn shift_octaves(spectrum: &mut [Complex32], u: &Unitary, pv: &mut Vocoder, mix: f32, span: usize) {
    // Phase a bin-centred partial advances per hop, per bin of frequency
    let hop_phase = TAU * HOP as f32 / FFT_SIZE as f32;
    pv.moved.fill(Complex32::new(0.0, 0.0));
    pv.loudest.fill(0.0);
    pv.stay.fill(1.0);
    pv.next_phase.copy_from_slice(&pv.out_phase);

    pv.peaks.clear();
    for k in OCTAVE_BASE..=TOP_BIN {
        let m = spectrum[k].norm_sqr();
        let louder = |j: usize| spectrum[j].norm_sqr() > m;
        let left = (k.saturating_sub(2).max(OCTAVE_BASE)..k)
            .any(|j| louder(j) || spectrum[j].norm_sqr() == m);
        let right = (k + 1..=(k + 2).min(TOP_BIN)).any(louder);
        if m > 0.0 && !left && !right {
            pv.peaks.push(k);
        }
    }

    for (i, &p) in pv.peaks.iter().enumerate() {
        let start = if i == 0 {
            OCTAVE_BASE
        } else {
            (pv.peaks[i - 1] + p) / 2 + 1
        };
        let end = pv.peaks.get(i + 1).map_or(TOP_BIN, |&q| (p + q) / 2);
        let phase = spectrum[p].arg();
        let deviation = phase - pv.prev_phase[p] - hop_phase * p as f32;
        // The partial's true frequency in bins, from how far its phase drifted off the bin centre.
        // A sinusoid's loudest bin is within half a bin of it; noise peaks can claim more, and
        // near the lobe's edge the amplitude below would blow up.
        let offset = (deviation - TAU * (deviation / TAU).round()) / hop_phase;
        let freq = p as f32 + offset.clamp(-0.5, 0.5);
        // Half the partial's amplitude (the positive-frequency half of the sinusoid), capped so a
        // redrawn copy never holds more energy than the region it stands for. Clean partials
        // aren't affected; noise is, and this keeps the moved energy within the Born-rule share.
        let model = spectrum[p].norm_sqr() / window_at(&pv.window, p as f32 - freq).norm_sqr();
        let region: f32 = spectrum[start..=end].iter().map(|x| x.norm_sqr()).sum();
        let amp = (model.min(region / WINDOW_ENERGY)).sqrt();
        let c = octave(p);
        let reach = |r: usize| r.abs_diff(c) <= span;
        let kept: f32 = (0..u.dim)
            .filter(|&r| reach(r))
            .map(|r| u.cols[c][r].norm_sqr())
            .sum();
        // kept > 0 unless the circuit sends this octave entirely out of reach; then leave it be
        let scale = if kept > 1e-12 {
            kept.sqrt().recip()
        } else {
            0.0
        };
        for r in (0..u.dim).filter(|&r| r != c && reach(r)) {
            let gain = u.cols[c][r].norm() * scale;
            let to = freq * 2f32.powi(r as i32 - c as i32);
            let target = to.round() as usize;
            if gain == 0.0 || !(1..=TOP_BIN).contains(&target) {
                continue;
            }
            let out = (pv.out_phase[target] + hop_phase * to).rem_euclid(TAU);
            if gain * amp > pv.loudest[target] {
                pv.loudest[target] = gain * amp;
                pv.next_phase[target] = out;
            }
            let partial = Complex32::from_polar(gain * amp, out);
            for k in target.saturating_sub(LOBE).max(1)..=(target + LOBE).min(TOP_BIN) {
                pv.moved[k] += partial * window_at(&pv.window, k as f32 - to);
            }
        }
        let stay = if scale > 0.0 {
            u.cols[c][c].norm() * scale
        } else {
            1.0
        };
        pv.stay[start..=end].fill(stay);
    }

    // Copies landing on each other can add up in phase (noise moved down an octave piles into half
    // as many bins), so cap the frame at its input energy: Octaves mode must never get louder
    let (mut dry_energy, mut wet_energy) = (0.0, 0.0);
    for ((x, stay), moved) in spectrum
        .iter()
        .zip(&pv.stay)
        .zip(&pv.moved)
        .take(TOP_BIN + 1)
        .skip(1)
    {
        dry_energy += x.norm_sqr();
        wet_energy += (x * stay + moved).norm_sqr();
    }
    let cap = if wet_energy > dry_energy {
        (dry_energy / wet_energy).sqrt()
    } else {
        1.0
    };
    for (k, bin) in spectrum.iter_mut().enumerate().take(TOP_BIN + 1).skip(1) {
        let x = *bin;
        pv.prev_phase[k] = x.arg();
        let wet = (x * pv.stay[k] + pv.moved[k]) * cap;
        *bin = x + (wet - x) * mix;
    }
    std::mem::swap(&mut pv.out_phase, &mut pv.next_phase);
}

/// Finds attacks from jumps in high-frequency content (bin-weighted energy) between frames.
#[derive(Clone, Copy, Default)]
struct Detector {
    prev_hfc: f32,
    /// 0..1, jumps up on an attack and decays
    env: f32,
}

impl Detector {
    fn next(&mut self, spectrum: &[Complex32], decay: f32) -> f32 {
        let hfc: f32 = spectrum
            .iter()
            .enumerate()
            .map(|(k, x)| k as f32 * x.norm_sqr())
            .sum();
        let rise_db = 10.0 * (hfc.max(HFC_FLOOR) / self.prev_hfc.max(HFC_FLOOR)).log10();
        self.prev_hfc = hfc;
        let onset = ((rise_db - ONSET_DB) / ONSET_RANGE_DB).clamp(0.0, 1.0);
        self.env = (self.env * decay).max(onset);
        self.env
    }
}

struct Tangle {
    params: Arc<TangleParams>,
    sample_rate: f32,
    stft: util::StftHelper,
    /// sqrt-Hann, used for both analysis and synthesis
    window: Vec<f32>,
    r2c: Arc<dyn RealToComplex<f32>>,
    c2r: Arc<dyn ComplexToReal<f32>>,
    spectrum: Vec<Complex32>,
    /// Left (and mono), and the circuit Spread morphs the right channel toward
    base_angles: [Angles; 2],
    unitaries: [Unitary; 2],
    /// In cycles of the drift rate
    drift_phase: f64,
    detectors: [Detector; 2],
    vocoders: [Vocoder; 2],
}

#[derive(Params)]
struct TangleParams {
    #[persist = "editor-size"]
    editor_state: Arc<EguiState>,

    #[id = "mode"]
    mode: EnumParam<Mode>,
    #[id = "amount"]
    amount: FloatParam,
    #[id = "mix"]
    mix: FloatParam,
    #[id = "width"]
    width: FloatParam,
    #[id = "qubits"]
    qubits: IntParam,
    #[id = "layers"]
    layers: IntParam,
    #[id = "twist"]
    twist: FloatParam,
    #[id = "coupling"]
    coupling: FloatParam,
    #[id = "drift"]
    drift: FloatParam,
    #[id = "drift-rate"]
    drift_rate: FloatParam,
    #[id = "transient"]
    transient: FloatParam,
    #[id = "spread"]
    spread: FloatParam,
    #[id = "span"]
    span: IntParam,
}

fn knob(name: &str, default: f32) -> FloatParam {
    FloatParam::new(name, default, FloatRange::Linear { min: 0.0, max: 1.0 })
        // Read once per hop (about 11 ms), so this ramps over a handful of frames
        .with_smoother(SmoothingStyle::Linear(50.0))
        .with_unit("%")
        .with_value_to_string(formatters::v2s_f32_percentage(0))
        .with_string_to_value(formatters::s2v_f32_percentage())
}

impl Default for TangleParams {
    fn default() -> Self {
        Self {
            editor_state: EguiState::from_size(360, 330),
            mode: EnumParam::new("Mode", Mode::Bands),
            amount: knob("Amount", 0.5),
            mix: knob("Mix", 1.0),
            // Steps between band sizes, so no smoother: each change jumps (crossfaded by the overlap)
            width: knob("Width", 0.2).with_smoother(SmoothingStyle::None),
            qubits: IntParam::new(
                "Qubits",
                4,
                IntRange::Linear {
                    min: 3,
                    max: MAX_QUBITS as i32,
                },
            )
            .with_value_to_string(Arc::new(|n| format!("{n} ({} bands)", 1 << n)))
            .with_string_to_value(Arc::new(|s| s.split_whitespace().next()?.parse().ok())),
            layers: IntParam::new(
                "Layers",
                2,
                IntRange::Linear {
                    min: 1,
                    max: MAX_LAYERS as i32,
                },
            ),
            twist: FloatParam::new(
                "Twist",
                0.0,
                FloatRange::Linear {
                    min: -1.0,
                    max: 1.0,
                },
            )
            .with_smoother(SmoothingStyle::Linear(50.0))
            .with_unit("%")
            .with_value_to_string(formatters::v2s_f32_percentage(0))
            .with_string_to_value(formatters::s2v_f32_percentage()),
            coupling: knob("Coupling", 1.0),
            drift: knob("Drift", 0.25),
            drift_rate: FloatParam::new(
                "Drift Rate",
                0.2,
                FloatRange::Skewed {
                    min: 0.01,
                    max: 2.0,
                    factor: FloatRange::skew_factor(-2.0),
                },
            )
            .with_unit(" Hz")
            .with_value_to_string(formatters::v2s_f32_rounded(2)),
            transient: knob("Transient Preserve", 0.5),
            // 0 so projects saved before Spread existed sound the same
            spread: knob("Spread", 0.0),
            // 7 reaches every octave, the 0.3.0 sound
            span: IntParam::new("Octave Span", 7, IntRange::Linear { min: 1, max: 7 })
                .with_value_to_string(Arc::new(|n| format!("±{n} oct")))
                .with_string_to_value(Arc::new(|s| {
                    s.trim_start_matches('±')
                        .split_whitespace()
                        .next()?
                        .parse()
                        .ok()
                })),
        }
    }
}

impl Default for Tangle {
    fn default() -> Self {
        let mut planner = RealFftPlanner::new();
        let r2c = planner.plan_fft_forward(FFT_SIZE);
        Self {
            params: Arc::new(TangleParams::default()),
            sample_rate: 48000.0,
            stft: util::StftHelper::new(2, FFT_SIZE, 0),
            // Periodic, so the squared windows sum to a constant at any hop of FFT_SIZE / 2^k
            window: (0..FFT_SIZE)
                .map(|n| (PI * n as f32 / FFT_SIZE as f32).sin())
                .collect(),
            spectrum: r2c.make_output_vec(),
            c2r: planner.plan_fft_inverse(FFT_SIZE),
            r2c,
            base_angles: [Angles::random(SEED), Angles::random(SEED + 1)],
            unitaries: std::array::from_fn(|_| Unitary::new(4, 0, &Angles::default())),
            drift_phase: 0.0,
            detectors: Default::default(),
            // Both channels up front, so a mono layout needs nothing new in initialize()
            vocoders: [Vocoder::new(), Vocoder::new()],
        }
    }
}

impl Tangle {
    /// `next` is called once per STFT frame for that frame's knob values.
    fn run(&mut self, buffer: &mut Buffer, mut next: impl FnMut() -> Settings) {
        // The FFT round trip scales by FFT_SIZE, and Hann windows at this overlap sum to OVERLAP / 2
        let gain = 2.0 / (FFT_SIZE * OVERLAP) as f32;
        let hop_seconds = HOP as f32 / self.sample_rate;
        let decay = (-hop_seconds / TRANSIENT_RELEASE).exp();
        let (window, r2c, c2r, spectrum) = (&self.window, &self.r2c, &self.c2r, &mut self.spectrum);
        let (bases, unitaries) = (&self.base_angles, &mut self.unitaries);
        let (drift_phase, detectors) = (&mut self.drift_phase, &mut self.detectors);
        let vocoders = &mut self.vocoders;
        let mut s = Settings::default();
        self.stft
            .process_overlap_add(buffer, OVERLAP, |channel, frame| {
                // Channels come in order for each frame
                if channel == 0 {
                    s = next();
                    *drift_phase += (s.drift_rate * hop_seconds) as f64;
                }
                let base = blend(
                    &bases[0],
                    &bases[1],
                    if channel == 0 { 0.0 } else { s.spread },
                );
                let qubits = if s.octaves { OCTAVE_QUBITS } else { s.qubits };
                let unitary = &mut unitaries[channel];
                unitary.build(qubits, s.layers, &s.angles(&base, *drift_phase));
                for (x, w) in frame.iter_mut().zip(window) {
                    *x *= w;
                }
                // Neither can fail: the buffer sizes come from the plans, and DC and Nyquist
                // aren't mixed so they stay real
                r2c.process_with_scratch(frame, spectrum, &mut []).unwrap();
                let attack = detectors[channel].next(spectrum, decay);
                let mix = s.mix * (1.0 - s.transient * attack);
                if s.octaves {
                    shift_octaves(spectrum, unitary, &mut vocoders[channel], mix, s.span);
                } else {
                    let k_len = s.band_bins(unitary.dim);
                    let block = k_len * unitary.dim;
                    let mixed = TOP_BIN / block * block;
                    mix_bands(&mut spectrum[1..=mixed], unitary, k_len, mix);
                }
                c2r.process_with_scratch(spectrum, frame, &mut []).unwrap();
                for (x, w) in frame.iter_mut().zip(window) {
                    *x *= w * gain;
                }
            });
    }
}

impl Plugin for Tangle {
    const NAME: &'static str = "Tangle";
    const VENDOR: &'static str = "accraze";
    const URL: &'static str = "";
    const EMAIL: &'static str = "";
    const VERSION: &'static str = env!("CARGO_PKG_VERSION");

    const AUDIO_IO_LAYOUTS: &'static [AudioIOLayout] = &[
        AudioIOLayout {
            main_input_channels: NonZeroU32::new(2),
            main_output_channels: NonZeroU32::new(2),
            ..AudioIOLayout::const_default()
        },
        AudioIOLayout {
            main_input_channels: NonZeroU32::new(1),
            main_output_channels: NonZeroU32::new(1),
            ..AudioIOLayout::const_default()
        },
    ];

    const SAMPLE_ACCURATE_AUTOMATION: bool = true;

    type SysExMessage = ();
    type BackgroundTask = ();

    fn params(&self) -> Arc<dyn Params> {
        self.params.clone()
    }

    fn editor(&mut self, _async_executor: AsyncExecutor<Self>) -> Option<Box<dyn Editor>> {
        // shortcut: nih-plug's generic sliders, swap for a hand-painted editor.rs like Spool's once the controls settle
        let params = self.params.clone();
        create_egui_editor(
            self.params.editor_state.clone(),
            (),
            |_, _| {},
            move |ctx, setter, _| {
                egui::CentralPanel::default().show(ctx, |ui| {
                    generic_ui::create(ui, params.clone(), setter, GenericSlider);
                });
            },
        )
    }

    fn initialize(
        &mut self,
        audio_io_layout: &AudioIOLayout,
        buffer_config: &BufferConfig,
        context: &mut impl InitContext<Self>,
    ) -> bool {
        self.sample_rate = buffer_config.sample_rate;
        let channels = audio_io_layout
            .main_output_channels
            .map_or(0, NonZeroU32::get);
        // Allocate here: reset() can run on the audio thread
        self.stft = util::StftHelper::new(channels as usize, FFT_SIZE, 0);
        context.set_latency_samples(self.stft.latency_samples());
        true
    }

    fn reset(&mut self) {
        // Also clears the buffers
        self.stft.set_block_size(FFT_SIZE);
        self.drift_phase = 0.0;
        self.detectors = Default::default();
        for pv in &mut self.vocoders {
            pv.reset();
        }
    }

    fn process(
        &mut self,
        buffer: &mut Buffer,
        _aux: &mut AuxiliaryBuffers,
        _context: &mut impl ProcessContext<Self>,
    ) -> ProcessStatus {
        let p = self.params.clone();
        let step = HOP as u32;
        self.run(buffer, || Settings {
            amount: p.amount.smoothed.next_step(step),
            mix: p.mix.smoothed.next_step(step),
            octaves: p.mode.value() == Mode::Octaves,
            width: p.width.value(),
            qubits: p.qubits.value() as usize,
            layers: p.layers.value() as usize,
            twist: p.twist.smoothed.next_step(step),
            coupling: p.coupling.smoothed.next_step(step),
            drift: p.drift.smoothed.next_step(step),
            drift_rate: p.drift_rate.value(),
            transient: p.transient.smoothed.next_step(step),
            spread: p.spread.smoothed.next_step(step),
            span: p.span.value() as usize,
        });
        ProcessStatus::Normal
    }
}

impl ClapPlugin for Tangle {
    const CLAP_ID: &'static str = "com.accraze.tangle";
    const CLAP_DESCRIPTION: Option<&'static str> =
        Some("Spectral smear through a simulated quantum circuit");
    const CLAP_MANUAL_URL: Option<&'static str> = None;
    const CLAP_SUPPORT_URL: Option<&'static str> = None;
    const CLAP_FEATURES: &'static [ClapFeature] = &[
        ClapFeature::AudioEffect,
        ClapFeature::Stereo,
        ClapFeature::Mono,
        ClapFeature::Filter,
    ];
}

impl Vst3Plugin for Tangle {
    const VST3_CLASS_ID: [u8; 16] = *b"TangleSpectralFx";
    const VST3_SUBCATEGORIES: &'static [Vst3SubCategory] =
        &[Vst3SubCategory::Fx, Vst3SubCategory::Filter];
}

nih_export_clap!(Tangle);
nih_export_vst3!(Tangle);

#[cfg(test)]
mod tests {
    use super::*;

    /// Full Amount, every modifier off
    const FULL: Settings = Settings {
        amount: 1.0,
        mix: 1.0,
        width: 1.0,
        qubits: 4,
        layers: 2,
        twist: 0.0,
        coupling: 1.0,
        drift: 0.0,
        drift_rate: 0.2,
        transient: 0.0,
        spread: 0.0,
        octaves: false,
        span: 7,
    };
    const OCTAVES: Settings = Settings {
        octaves: true,
        ..FULL
    };

    /// Runs `input` through the plugin as a host would, in uneven blocks that don't line up
    /// with the hop, and returns the output.
    fn render(input: &[Vec<f32>], settings: Settings) -> Vec<Vec<f32>> {
        let mut plugin = Tangle::default();
        let mut out = input.to_vec();
        let mut start = 0;
        while start < out[0].len() {
            let end = (start + 333).min(out[0].len());
            let mut buffer = Buffer::default();
            // Safety: the slices borrow `out`, which outlives `buffer`
            unsafe {
                buffer.set_slices(end - start, |slices| {
                    slices.extend(out.iter_mut().map(|ch| &mut ch[start..end]))
                })
            };
            plugin.run(&mut buffer, || settings);
            start = end;
        }
        out
    }

    fn noise(len: usize, seed: u32) -> Vec<f32> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                (state >> 8) as f32 / (1 << 23) as f32 - 1.0
            })
            .collect()
    }

    fn energy(x: &[Complex32]) -> f32 {
        x.iter().map(|c| c.norm_sqr()).sum()
    }

    fn rms(x: &[f32]) -> f32 {
        (x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32).sqrt()
    }

    /// Largest difference between `input` and `output` delayed by the latency
    fn null_error(input: &[f32], output: &[f32]) -> f32 {
        input
            .iter()
            .zip(&output[FFT_SIZE..])
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f32::max)
    }

    #[test]
    fn zero_amount_or_mix_is_a_null() {
        let input = vec![noise(48_000, 1), noise(48_000, 2)];
        // Every other knob pushed, to check none of them leak through
        let busy = Settings {
            qubits: MAX_QUBITS,
            layers: MAX_LAYERS,
            width: 0.0,
            twist: 1.0,
            drift: 1.0,
            drift_rate: 2.0,
            transient: 1.0,
            spread: 1.0,
            span: 1,
            ..FULL
        };
        for s in [
            Settings {
                amount: 0.0,
                ..busy
            },
            Settings { mix: 0.0, ..busy },
            Settings {
                amount: 0.0,
                octaves: true,
                ..busy
            },
            Settings {
                mix: 0.0,
                octaves: true,
                ..busy
            },
        ] {
            let output = render(&input, s);
            for (inp, out) in input.iter().zip(&output) {
                assert!(out[..FFT_SIZE].iter().all(|&y| y.abs() < 1e-6));
                let err = null_error(inp, out);
                println!("null test max error: {err:e}");
                assert!(err < 1e-5, "max error {err}");
            }
        }
        assert_eq!(Tangle::default().stft.latency_samples() as usize, FFT_SIZE);
    }

    #[test]
    fn mixing_keeps_each_frames_energy() {
        let re = noise(TOP_BIN, 3);
        let im = noise(TOP_BIN, 4);
        let frame: Vec<Complex32> = re
            .iter()
            .zip(&im)
            .map(|(&r, &i)| Complex32::new(r, i))
            .collect();
        let base = Angles::random(SEED);
        for qubits in 3..=MAX_QUBITS {
            for layers in 1..=MAX_LAYERS {
                for width in [0.0, 0.5, 1.0] {
                    let s = Settings {
                        qubits,
                        layers,
                        width,
                        twist: 0.3,
                        drift: 0.5,
                        ..FULL
                    };
                    let u = Unitary::new(qubits, layers, &s.angles(&base, 1.7));
                    let k_len = s.band_bins(u.dim);
                    assert_eq!(k_len % OVERLAP, 0);
                    let mixed = TOP_BIN / (k_len * u.dim) * k_len * u.dim;
                    let mut bins = frame[..mixed].to_vec();
                    mix_bands(&mut bins, &u, k_len, 1.0);
                    let ratio = energy(&bins) / energy(&frame[..mixed]);
                    let case = format!("{qubits} qubits, {layers} layers, width {width}");
                    assert!((ratio - 1.0).abs() < 1e-4, "{case}: {ratio}");
                    assert!(bins != frame[..mixed], "{case} didn't mix");
                }
            }
        }
    }

    #[test]
    fn output_level_stays_close() {
        // Frames overlap, so per-frame energy preservation doesn't make the output level exact.
        // White noise should still come out at about the same RMS.
        let input = vec![noise(96_000, 5), noise(96_000, 6)];
        let cases = [
            ("full", FULL),
            ("narrow", Settings { width: 0.0, ..FULL }),
            (
                "fast drift",
                Settings {
                    drift: 1.0,
                    drift_rate: 2.0,
                    ..FULL
                },
            ),
        ];
        for (name, s) in cases {
            let output = render(&input, s);
            let db = 20.0 * (rms(&output[0][FFT_SIZE..]) / rms(&input[0])).log10();
            println!("output level, {name}: {db:+.2} dB");
            assert!(db.abs() < 1.0, "{name}: {db} dB");
        }
        // Octaves mode only moves what looks like partials, so noise loses some of its moved
        // share. It must never come out louder.
        for (name, s) in [
            ("octaves", OCTAVES),
            ("octaves, span 1", Settings { span: 1, ..OCTAVES }),
            (
                "octaves, 4 layers, twisted",
                Settings {
                    layers: 4,
                    twist: 0.7,
                    ..OCTAVES
                },
            ),
        ] {
            let output = render(&input, s);
            let db = 20.0 * (rms(&output[0][FFT_SIZE..]) / rms(&input[0])).log10();
            println!("output level on noise, {name}: {db:+.2} dB");
            assert!((-6.0..0.5).contains(&db), "{name}: {db} dB");
        }
    }

    #[test]
    fn twist_coupling_and_drift_change_the_circuit() {
        let base = Angles::random(SEED);
        let still = FULL.angles(&base, 0.0);
        let changed = [
            Settings { twist: 0.5, ..FULL }.angles(&base, 0.0),
            Settings {
                coupling: 0.5,
                ..FULL
            }
            .angles(&base, 0.0),
            Settings { drift: 0.5, ..FULL }.angles(&base, 0.0),
        ];
        for a in changed {
            assert_ne!(a, still);
        }
        let drifting = Settings { drift: 0.5, ..FULL };
        assert_ne!(drifting.angles(&base, 0.0), drifting.angles(&base, 0.3));
        assert_eq!(
            FULL.angles(&base, 0.0),
            FULL.angles(&base, 0.3),
            "no drift, no motion"
        );
    }

    #[test]
    fn transient_preserve_keeps_attacks_dry() {
        // Quiet noise, then a loud decaying burst like a pick attack
        let len = 48_000;
        let onset = 24_000;
        let burst = noise(len, 7);
        let quiet = noise(len, 8);
        let x: Vec<f32> = (0..len)
            .map(|n| {
                let pluck = if n >= onset {
                    burst[n] * (-((n - onset) as f32) / 4000.0).exp()
                } else {
                    0.0
                };
                0.001 * quiet[n] + 0.5 * pluck
            })
            .collect();
        let input = vec![x.clone(), x];
        // Error in the 2048 samples around the onset, against the dry (delayed) input
        let attack_error = |transient| {
            let out = render(&input, Settings { transient, ..FULL });
            let range = onset - 512..onset + 1536;
            let err: Vec<f32> = range.map(|n| input[0][n] - out[0][n + FFT_SIZE]).collect();
            rms(&err)
        };
        let (off, on) = (attack_error(0.0), attack_error(1.0));
        println!("attack error: {off:.4} without, {on:.4} with Transient Preserve");
        assert!(on < 0.5 * off, "{on} vs {off}");
    }

    #[test]
    fn spread_gives_each_channel_its_own_circuit() {
        let x = noise(48_000, 9);
        let input = vec![x.clone(), x];
        for base in [FULL, OCTAVES] {
            let same = render(&input, base);
            assert_eq!(
                same[0], same[1],
                "Spread 0 must keep the channels identical"
            );
            let spread = render(
                &input,
                Settings {
                    spread: 1.0,
                    ..base
                },
            );
            assert_eq!(
                spread[0], same[0],
                "Spread must leave the left channel alone"
            );
            let (l, r) = (&spread[0][FFT_SIZE..], &spread[1][FFT_SIZE..]);
            let diff: Vec<f32> = l.iter().zip(r).map(|(a, b)| a - b).collect();
            let side = rms(&diff) / rms(l);
            let db = 20.0 * (rms(r) / rms(&input[1])).log10();
            println!("spread 100%: side/left {side:.2}, right channel {db:+.2} dB");
            assert!(side > 0.5, "channels barely differ: {side}");
            // Same bounds as output_level_stays_close
            let ok = if base.octaves {
                (-6.0..0.5).contains(&db)
            } else {
                db.abs() < 1.0
            };
            assert!(ok, "right channel level {db} dB");
        }
    }

    /// Power of `x` at `freq` Hz (48 kHz), through a Hann window so nearby tones don't leak in
    fn tone_power(x: &[f32], freq: f32) -> f32 {
        let n = x.len() as f32;
        let sum: Complex32 = x
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                let w = 0.5 - 0.5 * (TAU * i as f32 / n).cos();
                Complex32::from_polar(v * w, -TAU * freq * i as f32 / 48_000.0)
            })
            .sum();
        sum.norm_sqr() / (n * n)
    }

    #[test]
    fn octaves_mode_moves_partials_by_exact_octaves() {
        let input: Vec<f32> = (0..96_000)
            .map(|n| 0.5 * (TAU * 440.0 * n as f32 / 48_000.0).sin())
            .collect();
        // Each octave gets the Born-rule share |U[r][c]|² of the 440 Hz partial (octave c = 2)
        let u = Unitary::new(
            OCTAVE_QUBITS,
            2,
            &OCTAVES.angles(&Angles::random(SEED), 0.0),
        );
        let dry = tone_power(&input[FFT_SIZE + 8192..], 440.0);
        // Octave Span 1 only keeps 220 and 880 Hz, rescaled so the shares still sum to 1
        for span in [7, 1] {
            let output = render(
                &[input.clone(), input.clone()],
                Settings { span, ..OCTAVES },
            );
            let out = &output[0][FFT_SIZE + 8192..];
            let kept: f32 = (1..=3).map(|r| u.cols[2][r].norm_sqr()).sum();
            for (r, f) in [110.0, 220.0, 440.0, 880.0, 1760.0, 3520.0, 7040.0, 14080.0]
                .into_iter()
                .enumerate()
            {
                let share = match span {
                    7 => u.cols[2][r].norm_sqr(),
                    _ if r.abs_diff(2) <= 1 => u.cols[2][r].norm_sqr() / kept,
                    _ => 0.0,
                };
                let (got, expect) = (tone_power(out, f), share * dry);
                println!("span {span}, {f} Hz: {got:.3e}, expected {expect:.3e}");
                assert!(
                    (got - expect).abs() < 0.02 * expect + 1e-6 * dry,
                    "span {span}, {f} Hz: {got} vs {expect}"
                );
            }
        }
        let output = render(&[input.clone(), input.clone()], OCTAVES);
        // Skip the latency and the vocoder settling in
        let out = &output[0][FFT_SIZE + 8192..];
        let sum = |freqs: &[f32]| freqs.iter().map(|&f| tone_power(out, f)).sum::<f32>();
        let octaves = [110.0, 220.0, 880.0, 1760.0, 3520.0, 7040.0, 14080.0];
        // A fifth off, and 35 cents off (where copies would land without the vocoder)
        let fifths = octaves.map(|f| f * 1.5);
        let detuned = octaves.map(|f| f * 2f32.powf(-35.0 / 1200.0));
        let (on, off, flat) = (sum(&octaves), sum(&fifths), sum(&detuned));
        println!("octave copies {on:.2e}, at fifths {off:.2e}, 35 cents flat {flat:.2e}");
        assert!(on > 1e-4, "nothing moved: {on}");
        assert!(on > 100.0 * off, "{on} vs fifths {off}");
        assert!(on > 10.0 * flat, "{on} vs detuned {flat}");
        // Tonal input keeps its level: the vocoder keeps moved partials phase-coherent. (Noise
        // has no stable partials, so its moved share adds incoherently and loses up to ~4 dB.)
        let rms = |x: &[f32]| (x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32).sqrt();
        let db = 20.0 * (rms(out) / (0.5 / 2f32.sqrt())).log10();
        println!("output level, octaves on a sine: {db:+.2} dB");
        assert!(db.abs() < 1.0, "{db} dB");
    }
}
