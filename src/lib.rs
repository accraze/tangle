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
//! 4. Inverse FFT and overlap-add.

use nih_plug::prelude::*;
use nih_plug_egui::widgets::generic_ui::{self, GenericSlider};
use nih_plug_egui::{create_egui_editor, egui, EguiState};
use realfft::num_complex::Complex32;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::f32::consts::PI;
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
    base_angles: Angles,
    unitary: Unitary,
    /// In cycles of the drift rate
    drift_phase: f64,
    detectors: [Detector; 2],
}

#[derive(Params)]
struct TangleParams {
    #[persist = "editor-size"]
    editor_state: Arc<EguiState>,

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
            amount: knob("Amount", 0.5),
            mix: knob("Mix", 1.0),
            // Steps between band sizes, so no smoother: each change jumps (crossfaded by the overlap)
            width: knob("Width", 1.0).with_smoother(SmoothingStyle::None),
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
            base_angles: Angles::random(SEED),
            unitary: Unitary::new(4, 0, &Angles::default()),
            drift_phase: 0.0,
            detectors: Default::default(),
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
        let (base, unitary) = (&self.base_angles, &mut self.unitary);
        let (drift_phase, detectors) = (&mut self.drift_phase, &mut self.detectors);
        let mut s = Settings::default();
        self.stft
            .process_overlap_add(buffer, OVERLAP, |channel, frame| {
                // Channels come in order for each frame, and all share one U
                if channel == 0 {
                    s = next();
                    *drift_phase += (s.drift_rate * hop_seconds) as f64;
                    unitary.build(s.qubits, s.layers, &s.angles(base, *drift_phase));
                }
                for (x, w) in frame.iter_mut().zip(window) {
                    *x *= w;
                }
                // Neither can fail: the buffer sizes come from the plans, and DC and Nyquist
                // aren't mixed so they stay real
                r2c.process_with_scratch(frame, spectrum, &mut []).unwrap();
                let attack = detectors[channel].next(spectrum, decay);
                let k_len = s.band_bins(unitary.dim);
                let block = k_len * unitary.dim;
                let mixed = TOP_BIN / block * block;
                let mix = s.mix * (1.0 - s.transient * attack);
                mix_bands(&mut spectrum[1..=mixed], unitary, k_len, mix);
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
            width: p.width.value(),
            qubits: p.qubits.value() as usize,
            layers: p.layers.value() as usize,
            twist: p.twist.smoothed.next_step(step),
            coupling: p.coupling.smoothed.next_step(step),
            drift: p.drift.smoothed.next_step(step),
            drift_rate: p.drift_rate.value(),
            transient: p.transient.smoothed.next_step(step),
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
            ..FULL
        };
        for s in [
            Settings {
                amount: 0.0,
                ..busy
            },
            Settings { mix: 0.0, ..busy },
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
}
