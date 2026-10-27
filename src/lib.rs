//! A spectral smear for acoustic instruments.
//!
//! Signal chain:
//! 1. STFT: `FFT_SIZE` frames every `FFT_SIZE / OVERLAP` samples, periodic sqrt-Hann windows on
//!    analysis and synthesis, so overlap-add reconstructs the input exactly. nih-plug's
//!    `StftHelper` does the buffering; its latency is one frame.
//! 2. Spectral processing (identity for now).
//! 3. Inverse FFT and overlap-add.

use nih_plug::prelude::*;
use realfft::num_complex::Complex32;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::f32::consts::PI;
use std::sync::Arc;

const FFT_SIZE: usize = 2048;
const OVERLAP: usize = 4;

struct Tangle {
    params: Arc<TangleParams>,
    stft: util::StftHelper,
    /// sqrt-Hann, used for both analysis and synthesis
    window: Vec<f32>,
    r2c: Arc<dyn RealToComplex<f32>>,
    c2r: Arc<dyn ComplexToReal<f32>>,
    spectrum: Vec<Complex32>,
}

#[derive(Params, Default)]
struct TangleParams {}

impl Default for Tangle {
    fn default() -> Self {
        let mut planner = RealFftPlanner::new();
        let r2c = planner.plan_fft_forward(FFT_SIZE);
        Self {
            params: Arc::default(),
            stft: util::StftHelper::new(2, FFT_SIZE, 0),
            // Periodic, so the squared windows sum to a constant at any hop of FFT_SIZE / 2^k
            window: (0..FFT_SIZE)
                .map(|n| (PI * n as f32 / FFT_SIZE as f32).sin())
                .collect(),
            spectrum: r2c.make_output_vec(),
            c2r: planner.plan_fft_inverse(FFT_SIZE),
            r2c,
        }
    }
}

impl Tangle {
    fn run(&mut self, buffer: &mut Buffer) {
        // The FFT round trip scales by FFT_SIZE, and Hann windows at this overlap sum to OVERLAP / 2
        let gain = 2.0 / (FFT_SIZE * OVERLAP) as f32;
        let (window, r2c, c2r, spectrum) = (&self.window, &self.r2c, &self.c2r, &mut self.spectrum);
        self.stft
            .process_overlap_add(buffer, OVERLAP, |_channel, frame| {
                for (x, w) in frame.iter_mut().zip(window) {
                    *x *= w;
                }
                // Neither can fail: the buffer sizes come from the plans, and the identity leaves
                // the DC and Nyquist bins real
                r2c.process_with_scratch(frame, spectrum, &mut []).unwrap();
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
        self.run(buffer);
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

    /// Runs `input` through the plugin as a host would, in uneven blocks that don't line up
    /// with the hop, and returns the output.
    fn render(plugin: &mut Tangle, input: &[Vec<f32>]) -> Vec<Vec<f32>> {
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
            plugin.run(&mut buffer);
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

    #[test]
    fn identity_is_a_null() {
        let mut plugin = Tangle::default();
        let latency = plugin.stft.latency_samples() as usize;
        assert_eq!(latency, FFT_SIZE);
        let input = vec![noise(48_000, 1), noise(48_000, 2)];
        let output = render(&mut plugin, &input);
        for (inp, out) in input.iter().zip(&output) {
            assert!(out[..latency].iter().all(|&y| y.abs() < 1e-6));
            let err = inp
                .iter()
                .zip(&out[latency..])
                .map(|(x, y)| (x - y).abs())
                .fold(0.0, f32::max);
            println!("null test max error: {err:e}");
            assert!(err < 1e-5, "max error {err}");
        }
    }
}
