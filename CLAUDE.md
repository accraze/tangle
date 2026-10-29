# Tangle

Spectral smear plugin (CLAP/VST3) in Rust using nih-plug. Sibling of `../rodent`, `../wailer` and `../spool`
(same toolchain, CI and release process), but a workspace: the circuit simulator is its own crate in `core/`,
everything else is in `src/lib.rs`. The user runs it in Bitwig Studio on Linux (CLAP).

## Commands

```sh
cargo xtask bundle tangle --release && cp -r target/bundled/Tangle.clap ~/.clap/   # build + install
cargo test --release --workspace --lib -- --nocapture                             # plugin + core tests
cargo clippy --release --workspace --lib
```

- The system has an old distro Rust (1.75) in `/usr/bin` that comes before `~/.cargo/bin` on PATH and
  can't build nih-plug. Prefix commands with `PATH="$HOME/.cargo/bin:$PATH"`.
- Bitwig needs the plugin re-added or a restart to pick up a new build.
- Add user-facing changes to `CHANGELOG.md` under `[Unreleased]` as you make them.
- Releases: follow `docs/RELEASING.md` (CI publishes on annotated `v*` tags; never move an existing tag).
- The "trailing semicolon in macro" warnings come from nih-plug's export macros, not our code.

## DSP structure

- `core/src/lib.rs`: `Angles` (ry, rz, ent per layer and qubit) and `Unitary::build`, which simulates the
  circuit on each basis state into fixed `MAX_DIM` arrays (no allocation). Entanglers are controlled-RX, not
  CNOTs, so all-zero angles are exactly the identity. `Angles::random(seed)` is splitmix64.
- `src/lib.rs`: nih-plug's `StftHelper` (FFT 2048, overlap 4) with periodic sqrt-Hann windows. nih-plug's
  `util::window::hann` is the symmetric kind and doesn't reconstruct exactly, so don't swap it in.
- `Settings` holds one frame's knob values. `process()` reads it once per hop (`smoothed.next_step(HOP)`), the
  tests build it directly. `Settings::angles` maps knobs to circuit angles: everything is multiplied by Amount,
  which keeps Amount 0 an exact null. Keep it that way for new angle modifiers.
- `SEED` sets the base angles: changing it changes the sound of saved projects.
- U is rebuilt on channel 0 of each frame and shared by all channels. `Detector` (Transient Preserve) is per
  channel and scales the wet blend, not U.
- Band size (`Settings::band_bins`) must stay a multiple of `OVERLAP`: otherwise moved bins lose phase coherence
  across hops and the output drops ~4 dB (`output_level_stays_close` catches it). Never mix DC or Nyquist (bin
  1024): they must stay real or the inverse FFT's `unwrap()` panics.
- Octaves mode (`shift_octaves`): 3 qubits, U's index = octave of bins `OCTAVE_BASE·2^b..`. Peaks are partials;
  |U[r][c]|² of each moves to `freq·2^(r-c)` and is redrawn from the window's spectrum (`window_response`, tabled
  in `Vocoder::window` because the closed form costs 5x the CPU). Pitfalls the tests caught, keep them fixed:
  shifting lobes by whole bins detunes copies; a target bin's `out_phase` must be owned by its loudest partial
  (noise-floor peaks otherwise detune loud ones); the frequency estimate is clamped to ±0.5 bin and each copy's
  energy capped at its region's, or noise comes out 14 dB louder. `octaves_mode_moves_partials_by_exact_octaves`
  checks every octave's share against |U|², and the noise level check must stay below 0 dB.
  Octave Span drops shares more than `span` octaves away and rescales the rest by 1/sqrt(kept) so
  each column still sums to 1. Each frame is also capped at its input energy (copies can pile up in phase).
- The tests drive `run()` through a hand-built `Buffer` with stereo input, because `StftHelper` is sized for 2
  channels by default (mono only comes from `initialize()`).
- The audio thread must not allocate (`assert_process_allocs` is enabled).
- `CLAP_ID` and `VST3_CLASS_ID` must never change, or saved Bitwig projects lose the plugin.
