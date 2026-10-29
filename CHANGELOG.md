# Changelog

All notable changes to Tangle are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.4.0] - 2026-10-28

### Added
- Octave Span (Octaves mode): limits how far a partial can move, from ±1 to ±7 octaves. The
  shares within reach are scaled up so the level holds. ±7, the default, sounds the same as 0.3.0.

## [0.3.0] - 2026-10-28

### Added
- Spread: morphs the right channel toward its own circuit, so a mono source comes out wide.
  Each side is still energy-preserving, and 0% (the default) sounds the same as 0.2.0.

### Fixed
- Octaves mode could come out up to about 1 dB louder on dense material with some circuits.
  Each frame is now capped at its input energy.

## [0.2.0] - 2026-10-27

### Added
- Mode: Octaves. The circuit tangles each partial with its own octaves (8 octaves from about
  94 Hz up) instead of mixing frequency bands: a harmonic shimmer rather than a ring modulator.
  A phase vocoder keeps the octave copies exactly in tune. Bands, the default, sounds the same as
  0.1.0.

## [0.1.0] - 2026-10-27

### Added
- Spectral smear: the spectrum is split into bands, and the bands are mixed through the unitary
  matrix of a small simulated quantum circuit. Energy-preserving, so it can't blow up the level,
  and exactly transparent at Amount 0.
- Amount, Mix, Width, Qubits (8, 16 or 32 bands) and Layers.
- Twist, Coupling and Drift (with Drift Rate) to reshape the circuit and keep it moving.
- Transient Preserve, which pulls the effect back on pick attacks.
- Reports 2048 samples of latency (43 ms at 48 kHz).
- Generic slider GUI.

[Unreleased]: https://github.com/accraze/tangle/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/accraze/tangle/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/accraze/tangle/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/accraze/tangle/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/accraze/tangle/releases/tag/v0.1.0
