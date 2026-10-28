# Tangle

[![CI](https://github.com/accraze/tangle/actions/workflows/ci.yml/badge.svg)](https://github.com/accraze/tangle/actions/workflows/ci.yml)

A spectral smear plugin (CLAP / VST3) built with [nih-plug](https://github.com/robbert-vdh/nih-plug),
made with fingerstyle guitar in mind. It splits the spectrum into bands and tangles them together
through a small simulated quantum circuit. Keep it low for a shimmering blur around the notes, or
push it for drifting, metallic, frequency-shifted textures.

| Control | Range | What it does |
|---------|-------|--------------|
| Mode    | Bands / Octaves | What gets tangled: blocks of frequency bands, or each partial's octaves (see below) |
| Amount  | 0–100% | How far the circuit's gates rotate. 0% is exactly transparent |
| Mix     | 0–100% | Dry/wet |
| Width   | 0–100% | (Bands) Band size. Low settings use narrow bands (about 100–200 Hz) in small groups repeated up the spectrum, so energy moves in fine steps within a few kHz: a blur. High settings use bands over a kHz wide: bigger, more metallic jumps. Default 20% |
| Qubits  | 3–5 | (Bands) Number of bands that get mixed together: 8, 16 or 32. Octaves mode always uses 3 |
| Layers  | 1–4 | Depth of the circuit. More layers tangle the bands more thoroughly |
| Twist   | ±100% | Offsets every gate angle: a "different circuit" knob |
| Coupling | 0–100% | Strength of the entangling gates. At 0% the qubits rotate independently |
| Drift   | 0–100% | Each gate angle wobbles on its own slow LFO, so the smear keeps moving |
| Drift Rate | 0.01–2 Hz | How fast Drift moves |
| Transient Preserve | 0–100% | Pulls the effect back on pick attacks, so the transient stays intact, then lets go over 100 ms |
| Spread  | 0–100% | Morphs the right channel toward a different circuit: a mono guitar comes out wide. 0% = both sides identical |

Tangle reports 2048 samples of latency (43 ms at 48 kHz), which hosts compensate for. Width,
Qubits and Layers switch between discrete settings, so automating them jumps rather than glides.
The GUI is nih-plug's generic sliders for now.

## How it works

1. **STFT.** 2048-point FFTs every 512 samples (75% overlap) with sqrt-Hann windows on both sides,
   so with the effect off, overlap-add rebuilds the input exactly.
2. **Bands.** The bins are grouped into blocks of 2^N contiguous bands of K bins each (N = Qubits,
   K set by Width). For each offset k within a band, the 2^N bins at that offset form a vector, one
   bin per band.
3. **The circuit.** Each vector is multiplied by a 2^N × 2^N matrix U. U comes from simulating an
   N-qubit circuit: in each layer, an RY and RZ rotation on every qubit, then a ring of
   controlled-RX gates (qubit q controls q+1). Running the circuit on each basis state gives U's
   columns. The rotation angles are fixed random values, offset by Twist and Drift, and all of
   them are scaled by Amount.
4. **Why a quantum circuit.** Any circuit's matrix is unitary, so it keeps each frame's energy:
   it moves energy between bands, with interference, but never adds or removes any. The effect
   can't blow up the level, and Amount 0 (all angles zero) is exactly the identity.
5. **Why K is a multiple of 4.** Moving a bin by a multiple of FFT size / hop (4) bins keeps its
   phase advancing consistently from frame to frame, so overlapping frames add up instead of
   partly cancelling. Each move is a clean frequency shift, and the output level holds.

What you hear: each partial gets copied, at various strengths and phases, to its counterparts in
the other bands, which are shifted in frequency by whole bands. With narrow bands that's a dense
cloud of copies in fine frequency steps around every note; Drift keeps it moving.

### Octaves mode

Here the 3 qubits index 8 octaves, from about 94 Hz up to the top of the spectrum (at 48 kHz),
and U moves energy between a partial and its octaves instead of between bands. That's harmonic, so
it sounds like a shimmer rather than a ring modulator.

- Each spectral peak is treated as a partial. A phase vocoder measures its exact frequency.
- |U[r][c]|² is the share of a partial's energy that moves from octave c to octave r: the Born
  rule, the probability of measuring the circuit in state r. Because U is unitary, the shares
  from each octave add up to exactly 1, so energy is still conserved.
- The share that stays keeps its original bins (so Amount 0 is still exactly transparent). The
  shares that move are redrawn at exactly 2, 4, ½... times the partial's frequency, with their
  phase carried from frame to frame, so the copies are in tune.
- Only peaks move, so noisy content (pick noise, string squeak) loses a little of its moved
  share instead of being shifted. It never comes out louder.

There's no quantum hardware or quantum speed-up involved. The circuit is a compact way to build
a family of unitary matrices with a handful of musical knobs.

## Build

Needs a recent Rust from [rustup](https://rustup.rs) (distro packages are usually too old).
The GUI (egui) also needs X11/GL headers:

```sh
sudo apt install libx11-xcb-dev libxcb1-dev libxcb-dri2-0-dev libxcb-icccm4-dev libxcursor-dev libxkbcommon-dev libxcb-shape0-dev libxcb-xfixes0-dev libgl-dev
```

```sh
cargo xtask bundle tangle --release
cargo test --workspace --lib
```

Bundles land in `target/bundled/`. The circuit simulator lives in its own crate, `core/`.

## Install (Linux)

Grab `tangle-<version>-linux.zip` from [Releases](https://github.com/accraze/tangle/releases) and unzip it, or build it yourself. Then:

```sh
cp target/bundled/Tangle.clap ~/.clap/
cp -r target/bundled/Tangle.vst3 ~/.vst3/   # optional
```

Then rescan plug-ins in your DAW.

## Releasing

Changes are listed in [CHANGELOG.md](CHANGELOG.md). See [docs/RELEASING.md](docs/RELEASING.md) for the process. In short: bump the version, push an annotated `v*` tag, and CI publishes the release.

## License

Tangle is licensed under the [GNU General Public License v3.0 or later](LICENSE).
nih-plug itself is ISC licensed, but its VST3 bindings are GPLv3, so VST3 builds must be GPL anyway.
