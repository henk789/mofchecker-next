//! EQeq charge-equilibration kernel.
//!
//! This is a faithful translation of the EQeq C++ implementation
//! (github.com/lsmo-epfl/EQeq, `src/main.cpp`: `GetJ`, `Qeq`, `RoundCharges`,
//! `DetermineReciprocalLatticeVectors`) into Rust, so that the equilibrated
//! charges reproduce the reference bit-for-bit at its 3-decimal output
//! precision. The upstream EQeq is GPLv2; this module and the `eqeq`
//! subpackage that calls it are therefore GPL-licensed (see
//! `py/mofchecker_next/eqeq/LICENSE`), unlike the MIT core of the project.
//!
//! The kernel is purely numeric: the Python caller supplies fractional
//! coordinates, the six cell parameters, and the per-atom electronegativity
//! `X` and hardness `J` (computed from the vendored ionization/charge-center
//! tables). The kernel builds the cell vectors with EQeq's convention,
//! assembles the equilibration matrix, solves it, and rounds the charges with
//! EQeq's residual-charge adjustment.

use std::f64::consts::PI;

// glibc's erfc, identical to the one the EQeq shared library was built against.
extern "C" {
    fn erfc(x: f64) -> f64;
}

// EQeq physical constant (verbatim from main.cpp). The hydrogen X/J special case
// and all element parameters are computed Python-side and passed in as arrays.
const K: f64 = 14.4; // 1/(4 pi eps0) in Angstrom * eV

/// Inputs for one charge-equilibration solve.
pub struct EqeqInput<'a> {
    /// Fractional coordinates, one [a, b, c] per atom.
    pub frac_coords: &'a [[f64; 3]],
    /// Cell lengths a, b, c (Angstrom).
    pub cell_lengths: [f64; 3],
    /// Cell angles alpha, beta, gamma (degrees).
    pub cell_angles: [f64; 3],
    /// Per-atom electronegativity X (eV), already charge-center shifted.
    pub electronegativity: &'a [f64],
    /// Per-atom hardness/idempotential J (eV).
    pub hardness: &'a [f64],
    /// Net total charge of the cell (EQeq uses 0).
    pub total_charge: f64,
    /// Coulomb scaling parameter lambda (dielectric screening).
    pub lambda: f64,
    /// Ewald splitting parameter eta.
    pub eta: f64,
    /// Real-space expansion cell count (mR; +/- per axis).
    pub real_cells: i32,
    /// Reciprocal-space expansion cell count (mK; +/- per axis).
    pub recip_cells: i32,
    /// Number of digits for the final charge rounding.
    pub charge_precision: i32,
}

struct Cell {
    av: [f64; 3],
    bv: [f64; 3],
    cv: [f64; 3],
    hv: [f64; 3],
    jv: [f64; 3],
    kv: [f64; 3],
    volume: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn mag(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Build the real- and reciprocal-space cell vectors using EQeq's convention
/// (a along x, b in the xy plane).
fn build_cell(lengths: [f64; 3], angles_deg: [f64; 3]) -> Cell {
    let [a_len, b_len, c_len] = lengths;
    let alpha = angles_deg[0] * PI / 180.0;
    let beta = angles_deg[1] * PI / 180.0;
    let gamma = angles_deg[2] * PI / 180.0;

    let av = [a_len, 0.0, 0.0];
    let bv = [b_len * gamma.cos(), b_len * gamma.sin(), 0.0];
    let cv0 = c_len * beta.cos();
    let cv1 = (c_len * b_len * alpha.cos() - bv[0] * cv0) / bv[1];
    let cv2 = (c_len * c_len - cv0 * cv0 - cv1 * cv1).sqrt();
    let cv = [cv0, cv1, cv2];

    // Reciprocal-lattice vectors (2*pi convention), as in EQeq.
    let mut crs = cross(bv, cv);
    let mut pf = 2.0 * PI / dot(av, crs);
    let hv = [pf * crs[0], pf * crs[1], pf * crs[2]];
    crs = cross(cv, av);
    pf = 2.0 * PI / dot(bv, crs);
    let jv = [pf * crs[0], pf * crs[1], pf * crs[2]];
    crs = cross(av, bv);
    pf = 2.0 * PI / dot(cv, crs);
    let kv = [pf * crs[0], pf * crs[1], pf * crs[2]];

    let volume = (av[0] * cross(bv, cv)[0]
        + av[1] * cross(bv, cv)[1]
        + av[2] * cross(bv, cv)[2])
        .abs();

    Cell { av, bv, cv, hv, jv, kv, volume }
}

/// Convert fractional to Cartesian coordinates using EQeq's convention.
fn to_cartesian(frac: &[[f64; 3]], cell: &Cell) -> Vec<[f64; 3]> {
    frac
        .iter()
        .map(|f| {
            [
                f[0] * cell.av[0] + f[1] * cell.bv[0] + f[2] * cell.cv[0],
                f[0] * cell.av[1] + f[1] * cell.bv[1] + f[2] * cell.cv[1],
                f[0] * cell.av[2] + f[1] * cell.bv[2] + f[2] * cell.cv[2],
            ]
        })
        .collect()
}

/// EQeq's orbital-overlap term for a separation `r`.
#[inline]
fn orbital_overlap(j_i: f64, j_j: f64, rab: f64, rab_sq: f64) -> f64 {
    let jij = (j_i * j_j).sqrt();
    let a = jij / K;
    (-(a * a * rab_sq)).exp() * (2.0 * a - a * a * rab - 1.0 / rab)
}

/// Translation of EQeq's `GetJ(i, j)` for the periodic Ewald method.
fn get_j(
    i: usize,
    j: usize,
    pos: &[[f64; 3]],
    hardness: &[f64],
    cell: &Cell,
    lambda: f64,
    eta: f64,
    m_r: i32,
    m_k: i32,
) -> f64 {
    if i == j {
        // Orbital energy term over images (exclude origin).
        let mut orbital = 0.0;
        for u in -m_r..=m_r {
            for v in -m_r..=m_r {
                for w in -m_r..=m_r {
                    if u == 0 && v == 0 && w == 0 {
                        continue;
                    }
                    let d = lattice_shift(u, v, w, cell);
                    let rab_sq = dot(d, d);
                    let rab = rab_sq.sqrt();
                    orbital += orbital_overlap(hardness[i], hardness[j], rab, rab_sq);
                }
            }
        }
        // Real-space Coulomb component (exclude origin).
        let mut alpha_star = 0.0;
        for u in -m_r..=m_r {
            for v in -m_r..=m_r {
                for w in -m_r..=m_r {
                    if u == 0 && v == 0 && w == 0 {
                        continue;
                    }
                    let d = lattice_shift(u, v, w, cell);
                    let rab = mag(d);
                    alpha_star += unsafe { erfc(rab / eta) } / rab;
                }
            }
        }
        // Reciprocal-space component (exclude origin).
        let mut beta_star = 0.0;
        for u in -m_k..=m_k {
            for v in -m_k..=m_k {
                for w in -m_k..=m_k {
                    if u == 0 && v == 0 && w == 0 {
                        continue;
                    }
                    let rlv = recip_shift(u, v, w, cell);
                    let h = mag(rlv);
                    let b = 0.5 * h * eta;
                    beta_star += 1.0 / (h * h) * (-b * b).exp();
                }
            }
        }
        beta_star *= 4.0 * PI / cell.volume;

        hardness[i]
            + lambda * (K / 2.0) * (alpha_star + beta_star + orbital - 2.0 / (eta * PI.sqrt()))
    } else {
        let dij = [
            pos[i][0] - pos[j][0],
            pos[i][1] - pos[j][1],
            pos[i][2] - pos[j][2],
        ];
        // Orbital energy term over all images (include origin).
        let mut orbital = 0.0;
        for u in -m_r..=m_r {
            for v in -m_r..=m_r {
                for w in -m_r..=m_r {
                    let s = lattice_shift(u, v, w, cell);
                    let d = [dij[0] + s[0], dij[1] + s[1], dij[2] + s[2]];
                    let rab_sq = dot(d, d);
                    let rab = rab_sq.sqrt();
                    orbital += orbital_overlap(hardness[i], hardness[j], rab, rab_sq);
                }
            }
        }
        // Real-space Coulomb component over all images (include origin).
        let mut alpha = 0.0;
        for u in -m_r..=m_r {
            for v in -m_r..=m_r {
                for w in -m_r..=m_r {
                    let s = lattice_shift(u, v, w, cell);
                    let d = [dij[0] + s[0], dij[1] + s[1], dij[2] + s[2]];
                    let rab = mag(d);
                    alpha += unsafe { erfc(rab / eta) } / rab;
                }
            }
        }
        // Reciprocal-space component (exclude origin), phase from the bare dij.
        let mut beta = 0.0;
        for u in -m_k..=m_k {
            for v in -m_k..=m_k {
                for w in -m_k..=m_k {
                    if u == 0 && v == 0 && w == 0 {
                        continue;
                    }
                    let rlv = recip_shift(u, v, w, cell);
                    let h = mag(rlv);
                    let b = 0.5 * h * eta;
                    beta += (rlv[0] * dij[0] + rlv[1] * dij[1] + rlv[2] * dij[2]).cos()
                        / (h * h)
                        * (-b * b).exp();
                }
            }
        }
        beta *= 4.0 * PI / cell.volume;

        lambda * (K / 2.0) * (alpha + beta + orbital)
    }
}

#[inline]
fn lattice_shift(u: i32, v: i32, w: i32, cell: &Cell) -> [f64; 3] {
    let (u, v, w) = (u as f64, v as f64, w as f64);
    [
        u * cell.av[0] + v * cell.bv[0] + w * cell.cv[0],
        u * cell.av[1] + v * cell.bv[1] + w * cell.cv[1],
        u * cell.av[2] + v * cell.bv[2] + w * cell.cv[2],
    ]
}

#[inline]
fn recip_shift(u: i32, v: i32, w: i32, cell: &Cell) -> [f64; 3] {
    let (u, v, w) = (u as f64, v as f64, w as f64);
    [
        u * cell.hv[0] + v * cell.jv[0] + w * cell.kv[0],
        u * cell.hv[1] + v * cell.jv[1] + w * cell.kv[1],
        u * cell.hv[2] + v * cell.jv[2] + w * cell.kv[2],
    ]
}

/// EQeq's round-half-away-from-zero.
#[inline]
fn round_half_away(num: f64) -> f64 {
    if num > 0.0 {
        (num + 0.5).floor()
    } else {
        (num - 0.5).ceil()
    }
}

/// Translation of EQeq's `RoundCharges`.
fn round_charges(q: &mut [f64], digits: i32) {
    let factor = 10f64.powi(digits);
    let mut qsum = 0.0;
    for value in q.iter_mut() {
        *value = round_half_away(*value * factor) / factor;
        qsum += *value;
    }
    if qsum == 0.0 {
        return;
    }
    let num_to_adjust = ((qsum * factor).abs() + 0.5) as i32;
    let sign = if qsum > 0.0 { -1.0 } else { 1.0 };
    for value in q.iter_mut().take(num_to_adjust.max(0) as usize) {
        *value += sign * (1.0 / factor);
    }
}

/// Equilibrate charges. Returns one charge per atom.
pub fn equilibrate(input: &EqeqInput) -> Result<Vec<f64>, String> {
    let n = input.frac_coords.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if input.electronegativity.len() != n || input.hardness.len() != n {
        return Err("X/J arrays must match atom count".to_string());
    }

    let cell = build_cell(input.cell_lengths, input.cell_angles);
    if !cell.volume.is_finite() || cell.volume <= 0.0 {
        return Err("cell has non-positive volume".to_string());
    }
    let pos = to_cartesian(input.frac_coords, &cell);

    // Precompute the full J-matrix once: A[i][j] = jmat[i-1][j] - jmat[i][j].
    let mut jmat = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            jmat[i * n + j] = get_j(
                i,
                j,
                &pos,
                input.hardness,
                &cell,
                input.lambda,
                input.eta,
                input.real_cells,
                input.recip_cells,
            );
        }
    }

    // Assemble A x = b exactly as EQeq's Qeq().
    let mut a = vec![0.0f64; n * n];
    let mut b = vec![0.0f64; n];
    for col in 0..n {
        a[col] = 1.0; // first row all ones
    }
    b[0] = input.total_charge;
    for i in 1..n {
        b[i] = input.electronegativity[i] - input.electronegativity[i - 1];
        for j in 0..n {
            a[i * n + j] = jmat[(i - 1) * n + j] - jmat[i * n + j];
        }
    }

    let mut q = solve_linear(a, b, n)?;
    round_charges(&mut q, input.charge_precision);
    Ok(q)
}

/// Solve a dense linear system by Gaussian elimination with partial pivoting.
/// The Qeq system is square with a unique solution, so this reproduces EQeq's
/// Householder solve to full precision.
fn solve_linear(mut a: Vec<f64>, mut b: Vec<f64>, dim: usize) -> Result<Vec<f64>, String> {
    if dim == 1 {
        if a[0].abs() < 1e-300 {
            return Err("singular 1x1 system".to_string());
        }
        return Ok(vec![b[0] / a[0]]);
    }
    for col in 0..dim {
        let mut pivot = col;
        let mut best = a[col * dim + col].abs();
        for row in (col + 1)..dim {
            let v = a[row * dim + col].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1e-300 {
            return Err("singular charge-equilibration matrix".to_string());
        }
        if pivot != col {
            for k in 0..dim {
                a.swap(col * dim + k, pivot * dim + k);
            }
            b.swap(col, pivot);
        }
        let diag = a[col * dim + col];
        for row in (col + 1)..dim {
            let factor = a[row * dim + col] / diag;
            if factor != 0.0 {
                for k in col..dim {
                    a[row * dim + k] -= factor * a[col * dim + k];
                }
                b[row] -= factor * b[col];
            }
        }
    }
    let mut x = vec![0.0f64; dim];
    for col in (0..dim).rev() {
        let mut acc = b[col];
        for k in (col + 1)..dim {
            acc -= a[col * dim + k] * x[k];
        }
        x[col] = acc / a[col * dim + col];
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic(frac: Vec<[f64; 3]>, x: Vec<f64>, j: Vec<f64>) -> Vec<f64> {
        let input = EqeqInput {
            frac_coords: &frac,
            cell_lengths: [15.0, 15.0, 15.0],
            cell_angles: [90.0, 90.0, 90.0],
            electronegativity: &x,
            hardness: &j,
            total_charge: 0.0,
            lambda: 1.2,
            eta: 50.0,
            real_cells: 2,
            recip_cells: 2,
            charge_precision: 3,
        };
        equilibrate(&input).unwrap()
    }

    #[test]
    fn single_atom_zero_charge() {
        let q = cubic(vec![[0.5, 0.5, 0.5]], vec![5.0], vec![10.0]);
        assert_eq!(q, vec![0.0]);
    }

    #[test]
    fn charges_sum_to_zero_after_rounding() {
        let q = cubic(
            vec![[0.1, 0.1, 0.1], [0.2, 0.1, 0.1], [0.1, 0.2, 0.1]],
            vec![13.6, 7.5, 7.5],
            vec![10.0, 12.0, 12.0],
        );
        let total: f64 = q.iter().sum();
        assert!(total.abs() < 1e-9, "rounded charges must sum to zero: {total}");
    }

    #[test]
    fn round_half_away_matches_eqeq() {
        assert_eq!(round_half_away(0.5), 1.0);
        assert_eq!(round_half_away(-0.5), -1.0);
        assert_eq!(round_half_away(2.4), 2.0);
        assert_eq!(round_half_away(-2.6), -3.0);
    }
}
