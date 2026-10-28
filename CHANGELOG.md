# Changelog

All notable changes to Tangle are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

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

[Unreleased]: https://github.com/accraze/tangle/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/accraze/tangle/releases/tag/v0.1.0
