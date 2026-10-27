//! A spectral smear for acoustic instruments.
//!
//! Signal chain:
//! 1. STFT: `FFT_SIZE` frames every `HOP` samples, periodic sqrt-Hann windows on analysis and
//!    synthesis, so overlap-add reconstructs the input exactly. nih-plug's `StftHelper` does the
//!    buffering; its latency is one frame.
//! 2. Once per frame, the knobs are read and the circuit's unitary U is rebuilt (`tangle_core`).
//!    Bins from 1 up to `TOP_BIN` are split into 2^qubits contiguous bands of K bins (`band_bins`;
//!    the few left over at the very top stay put), and for each offset
//!    k the vector of bin k from every band is multiplied by U. U is unitary, so this keeps the
//!    frame's energy; Mix blends towards the identity in the spectrum, which equals a dry blend.
//! 3. Inverse FFT and overlap-add.

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

/// Knob values for one STFT frame
#[derive(Clone, Copy)]
struct Settings {
    amount: f32,
    mix: f32,
    qubits: usize,
    layers: usize,
}

impl Settings {
    /// Base angles scaled by Amount, so Amount 0 is exactly the identity.
    fn angles(&self, base: &Angles) -> Angles {
        let scale =
            |a: [[f32; MAX_QUBITS]; MAX_LAYERS]| a.map(|layer| layer.map(|x| x * self.amount));
        Angles {
            ry: scale(base.ry),
            rz: scale(base.rz),
            ent: scale(base.ent),
        }
    }
}

/// Bins per band. A multiple of OVERLAP, so moving a bin by whole bands shifts it by a multiple
/// of FFT_SIZE / HOP bins: its phase still advances consistently from hop to hop, overlapping frames
/// add up coherently, and the move is a clean frequency shift. Otherwise the frames partly cancel
/// (about -4 dB at full Amount) and the result gets phasey.
fn band_bins(dim: usize) -> usize {
    TOP_BIN / dim / OVERLAP * OVERLAP
}

/// Multiplies each cross-band vector of `bins` by `u`, blended with the input by `mix`.
fn mix_bands(bins: &mut [Complex32], u: &Unitary, mix: f32) {
    let k_len = bins.len() / u.dim;
    for k in 0..k_len {
        let mut x = [Complex32::new(0.0, 0.0); MAX_DIM];
        for (b, x) in x.iter_mut().enumerate().take(u.dim) {
            *x = bins[b * k_len + k];
        }
        for r in 0..u.dim {
            let ux: Complex32 = (0..u.dim).map(|c| u.cols[c][r] * x[c]).sum();
            bins[r * k_len + k] = x[r] + (ux - x[r]) * mix;
        }
    }
}

struct Tangle {
    params: Arc<TangleParams>,
    stft: util::StftHelper,
    /// sqrt-Hann, used for both analysis and synthesis
    window: Vec<f32>,
    r2c: Arc<dyn RealToComplex<f32>>,
    c2r: Arc<dyn ComplexToReal<f32>>,
    spectrum: Vec<Complex32>,
    base_angles: Angles,
    unitary: Unitary,
}

#[derive(Params)]
struct TangleParams {
    #[persist = "editor-size"]
    editor_state: Arc<EguiState>,

    #[id = "amount"]
    amount: FloatParam,
    #[id = "mix"]
    mix: FloatParam,
    #[id = "qubits"]
    qubits: IntParam,
    #[id = "layers"]
    layers: IntParam,
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
            editor_state: EguiState::from_size(360, 180),
            amount: knob("Amount", 0.5),
            mix: knob("Mix", 1.0),
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
        }
    }
}

impl Default for Tangle {
    fn default() -> Self {
        let mut planner = RealFftPlanner::new();
        let r2c = planner.plan_fft_forward(FFT_SIZE);
        Self {
            params: Arc::new(TangleParams::default()),
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
        }
    }
}

impl Tangle {
    /// `next` is called once per STFT frame for that frame's knob values.
    fn run(&mut self, buffer: &mut Buffer, mut next: impl FnMut() -> Settings) {
        // The FFT round trip scales by FFT_SIZE, and Hann windows at this overlap sum to OVERLAP / 2
        let gain = 2.0 / (FFT_SIZE * OVERLAP) as f32;
        let (window, r2c, c2r, spectrum) = (&self.window, &self.r2c, &self.c2r, &mut self.spectrum);
        let (base, unitary) = (&self.base_angles, &mut self.unitary);
        let mut mix = 0.0;
        self.stft
            .process_overlap_add(buffer, OVERLAP, |channel, frame| {
                // Channels come in order for each frame, and all share one U
                if channel == 0 {
                    let s = next();
                    unitary.build(s.qubits, s.layers, &s.angles(base));
                    mix = s.mix;
                }
                for (x, w) in frame.iter_mut().zip(window) {
                    *x *= w;
                }
                // Neither can fail: the buffer sizes come from the plans, and DC and Nyquist
                // aren't mixed so they stay real
                r2c.process_with_scratch(frame, spectrum, &mut []).unwrap();
                let mixed = band_bins(unitary.dim) * unitary.dim;
                mix_bands(&mut spectrum[1..=mixed], unitary, mix);
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
        _buffer_config: &BufferConfig,
        context: &mut impl InitContext<Self>,
    ) -> bool {
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
    }

    fn process(
        &mut self,
        buffer: &mut Buffer,
        _aux: &mut AuxiliaryBuffers,
        _context: &mut impl ProcessContext<Self>,
    ) -> ProcessStatus {
        let p = self.params.clone();
        self.run(buffer, || Settings {
            amount: p.amount.smoothed.next_step(HOP as u32),
            mix: p.mix.smoothed.next_step(HOP as u32),
            qubits: p.qubits.value() as usize,
            layers: p.layers.value() as usize,
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

    const FULL: Settings = Settings {
        amount: 1.0,
        mix: 1.0,
        qubits: 4,
        layers: 2,
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

    #[test]
    fn zero_amount_or_mix_is_a_null() {
        let input = vec![noise(48_000, 1), noise(48_000, 2)];
        let biggest = Settings {
            qubits: MAX_QUBITS,
            layers: MAX_LAYERS,
            ..FULL
        };
        let settings = [
            Settings {
                amount: 0.0,
                ..biggest
            },
            Settings {
                mix: 0.0,
                ..biggest
            },
        ];
        for s in settings {
            let output = render(&input, s);
            for (inp, out) in input.iter().zip(&output) {
                assert!(out[..FFT_SIZE].iter().all(|&y| y.abs() < 1e-6));
                let err = inp
                    .iter()
                    .zip(&out[FFT_SIZE..])
                    .map(|(x, y)| (x - y).abs())
                    .fold(0.0, f32::max);
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
                let s = Settings {
                    qubits,
                    layers,
                    ..FULL
                };
                let u = Unitary::new(qubits, layers, &s.angles(&base));
                let mut bins = frame.clone();
                mix_bands(&mut bins, &u, 1.0);
                let ratio = energy(&bins) / energy(&frame);
                assert!(
                    (ratio - 1.0).abs() < 1e-4,
                    "{qubits} qubits, {layers} layers: {ratio}"
                );
                assert!(bins != frame, "{qubits} qubits, {layers} layers didn't mix");
            }
        }
    }

    #[test]
    fn output_level_stays_close() {
        // Frames overlap, so per-frame energy preservation doesn't make the output level exact.
        // White noise should still come out at about the same RMS.
        let input = vec![noise(96_000, 5), noise(96_000, 6)];
        let output = render(&input, FULL);
        let rms = |x: &[f32]| (x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32).sqrt();
        let db = 20.0 * (rms(&output[0][FFT_SIZE..]) / rms(&input[0])).log10();
        println!("output level at full Amount: {db:+.2} dB");
        assert!(db.abs() < 1.0, "{db} dB");
    }
}
