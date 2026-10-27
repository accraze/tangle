//! A tiny statevector simulator that turns a layered quantum circuit into a unitary matrix.
//!
//! Each layer applies RY then RZ to every qubit, then a ring of controlled-RX gates (qubit q
//! controls q+1, wrapping). A controlled-RX at π is a CNOT up to phase; at 0 it's the identity,
//! so all-zero angles give exactly the identity matrix (a fixed CNOT ring would permute instead).
//!
//! Everything is sized for `MAX_QUBITS`, so building a unitary never allocates.

pub use num_complex::Complex32;
use std::f32::consts::PI;

pub const MAX_QUBITS: usize = 5;
pub const MAX_DIM: usize = 1 << MAX_QUBITS;
pub const MAX_LAYERS: usize = 4;

const ZERO: Complex32 = Complex32::new(0.0, 0.0);

/// Gate angles in radians, indexed `[layer][qubit]`. Default is all zero (the identity).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Angles {
    pub ry: [[f32; MAX_QUBITS]; MAX_LAYERS],
    pub rz: [[f32; MAX_QUBITS]; MAX_LAYERS],
    /// Controlled-RX from qubit q to q+1
    pub ent: [[f32; MAX_QUBITS]; MAX_LAYERS],
}

impl Angles {
    /// Rotations uniform in [-π, π) from `seed`, entanglers at π (full CNOTs).
    pub fn random(seed: u64) -> Self {
        let mut state = seed;
        let mut next = || {
            // splitmix64
            state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
            ((z >> 40) as f32 / (1u64 << 24) as f32 * 2.0 - 1.0) * PI
        };
        let mut a = Self::default();
        for l in 0..MAX_LAYERS {
            for q in 0..MAX_QUBITS {
                a.ry[l][q] = next();
                a.rz[l][q] = next();
                a.ent[l][q] = PI;
            }
        }
        a
    }
}

/// A `dim` x `dim` unitary, stored by column: `cols[c]` is the circuit applied to basis state c.
#[derive(Clone, Debug)]
pub struct Unitary {
    pub dim: usize,
    pub cols: [[Complex32; MAX_DIM]; MAX_DIM],
}

impl Unitary {
    pub fn new(qubits: usize, layers: usize, angles: &Angles) -> Self {
        let mut u = Self {
            dim: 0,
            cols: [[ZERO; MAX_DIM]; MAX_DIM],
        };
        u.build(qubits, layers, angles);
        u
    }

    /// Rebuilds in place, so the audio thread can call it once per hop.
    pub fn build(&mut self, qubits: usize, layers: usize, angles: &Angles) {
        assert!((1..=MAX_QUBITS).contains(&qubits) && layers <= MAX_LAYERS);
        self.dim = 1 << qubits;
        for (c, col) in self.cols.iter_mut().enumerate().take(self.dim) {
            *col = [ZERO; MAX_DIM];
            col[c] = Complex32::new(1.0, 0.0);
            run(&mut col[..self.dim], qubits, layers, angles);
        }
    }
}

fn run(state: &mut [Complex32], qubits: usize, layers: usize, a: &Angles) {
    for l in 0..layers {
        for q in 0..qubits {
            let (s, c) = (a.ry[l][q] / 2.0).sin_cos();
            apply_1q(
                state,
                q,
                None,
                [[c.into(), (-s).into()], [s.into(), c.into()]],
            );
            let (s, c) = (a.rz[l][q] / 2.0).sin_cos();
            let (lo, hi) = (Complex32::new(c, -s), Complex32::new(c, s));
            apply_1q(state, q, None, [[lo, ZERO], [ZERO, hi]]);
        }
        for q in 0..qubits {
            let (s, c) = (a.ent[l][q] / 2.0).sin_cos();
            let (c, s) = (Complex32::new(c, 0.0), Complex32::new(0.0, -s));
            apply_1q(state, (q + 1) % qubits, Some(q), [[c, s], [s, c]]);
        }
    }
}

/// Applies the 2x2 gate `g` to qubit `target`, only where qubit `control` (if any) is 1.
fn apply_1q(
    state: &mut [Complex32],
    target: usize,
    control: Option<usize>,
    g: [[Complex32; 2]; 2],
) {
    let t = 1 << target;
    for i in 0..state.len() {
        if i & t != 0 || control.is_some_and(|c| i & (1 << c) == 0) {
            continue;
        }
        let (a, b) = (state[i], state[i | t]);
        state[i] = g[0][0] * a + g[0][1] * b;
        state[i | t] = g[1][0] * a + g[1][1] * b;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_shapes() -> impl Iterator<Item = (usize, usize)> {
        (3..=MAX_QUBITS).flat_map(|n| (1..=MAX_LAYERS).map(move |l| (n, l)))
    }

    #[test]
    fn is_unitary() {
        for seed in 0..8 {
            for (n, layers) in all_shapes() {
                let u = Unitary::new(n, layers, &Angles::random(seed));
                for i in 0..u.dim {
                    for j in 0..u.dim {
                        // (U†U)_ij = <col i, col j>
                        let dot: Complex32 =
                            (0..u.dim).map(|r| u.cols[i][r].conj() * u.cols[j][r]).sum();
                        let expect = if i == j { 1.0 } else { 0.0 };
                        assert!(
                            (dot - expect).norm() < 1e-5,
                            "seed {seed}, {n} qubits, {layers} layers: (U†U)[{i}][{j}] = {dot}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn zero_angles_are_exactly_the_identity() {
        for (n, layers) in all_shapes() {
            let u = Unitary::new(n, layers, &Angles::default());
            for c in 0..u.dim {
                for r in 0..u.dim {
                    let expect = if r == c { 1.0 } else { 0.0 };
                    assert_eq!(
                        u.cols[c][r],
                        Complex32::new(expect, 0.0),
                        "{n} qubits, {layers} layers"
                    );
                }
            }
        }
    }

    #[test]
    fn same_seed_same_matrix() {
        let (a, b) = (
            Unitary::new(4, 2, &Angles::random(42)),
            Unitary::new(4, 2, &Angles::random(42)),
        );
        assert_eq!(a.cols, b.cols);
        assert_ne!(a.cols, Unitary::new(4, 2, &Angles::random(43)).cols);
    }

    #[test]
    fn entangler_ring_is_wired_in_order() {
        // No rotations, full entanglers: |q0=1> flips q1, which flips q2, which flips q0 back.
        // Each controlled-RX(π) flip adds a -i phase: |001> -> (-i)^3 |110> = i |110>.
        let mut a = Angles::default();
        a.ent[0] = [PI; MAX_QUBITS];
        let u = Unitary::new(3, 1, &a);
        assert!((u.cols[0b001][0b110] - Complex32::new(0.0, 1.0)).norm() < 1e-6);
        assert!(
            u.cols[0][0] == Complex32::new(1.0, 0.0),
            "|000> has no control set"
        );
    }
}
