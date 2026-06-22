// MIT License
//
// Copyright (c) 2026 Libre Ephemeris Contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

pub mod iau2000a;
pub mod iau2000b;

use crate::constants;
use crate::types::LeMat3;

/// Nutation components: longitude nutation (Δψ) and obliquity nutation (Δε) in radians.
#[derive(Debug, Clone, Copy)]
pub struct Nutation {
    pub dpsi: f64,
    pub deps: f64,
    pub eps0: f64, // mean obliquity
    pub eps: f64,  // true obliquity = eps0 + deps
    pub t_jcent: f64,
}

/// Compute nutation for a given Julian day and model ID.
pub fn compute_nutation(jd: f64, model: i32) -> Nutation {
    let t = (jd - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let eps0 = crate::transform::mean_obliquity_iau2006(t);

    let (dpsi, deps) = match model {
        constants::LE_NUT_IAU_2000A => iau2000a::nutation_lon_obl(t),
        constants::LE_NUT_IAU_2000B => iau2000b::nutation_lon_obl(t),
        constants::LE_NUT_IAU_1980 => iau2000b::nutation_lon_obl(t), // approximation
        constants::LE_NUT_IAU_1980_HERRING => iau2000b::nutation_lon_obl(t),
        constants::LE_NUT_WOOLARD_1953 => iau2000b::nutation_lon_obl(t),
        _ => iau2000b::nutation_lon_obl(t),
    };

    Nutation {
        dpsi,
        deps,
        eps0,
        eps: eps0 + deps,
        t_jcent: t,
    }
}

/// Build the nutation rotation matrix.
/// Accounts for both nutation in longitude and obliquity.
pub fn nutation_matrix(nut: &Nutation) -> LeMat3 {
    // Rotation sequence: R1(-eps0 - deps) * R3(-dpsi) * R1(eps0)
    let cos_eps0 = nut.eps0.cos();
    let sin_eps0 = nut.eps0.sin();
    let cos_eps = nut.eps.cos();
    let sin_eps = nut.eps.sin();
    let cos_dpsi = nut.dpsi.cos();
    let sin_dpsi = nut.dpsi.sin();

    LeMat3([
        [cos_dpsi, cos_eps0 * sin_dpsi, sin_eps0 * sin_dpsi],
        [-cos_eps * sin_dpsi, cos_eps * cos_eps0 * cos_dpsi + sin_eps * sin_eps0,
            cos_eps * sin_eps0 * cos_dpsi - sin_eps * cos_eps0],
        [-sin_eps * sin_dpsi, sin_eps * cos_eps0 * cos_dpsi - cos_eps * sin_eps0,
            sin_eps * sin_eps0 * cos_dpsi + cos_eps * cos_eps0],
    ])
}

/// Get combined precession-nutation matrix from J2000 to true equator of date.
pub fn precession_nutation_matrix(jd: f64, prec_model: i32, nut_model: i32) -> LeMat3 {
    let prec = crate::precession::precession_matrix_for_model(jd, prec_model);
    let nut = compute_nutation(jd, nut_model);
    let nmat = nutation_matrix(&nut);
    nmat.mul(&prec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nutation_j2000_reasonable() {
        let nut = compute_nutation(constants::LE_J2000, constants::LE_NUT_IAU_2000B);
        // Nutation at J2000 is ~14-18 arcseconds in longitude (~7-9e-5 rad)
        // and ~3-9 arcseconds in obliquity (~1.5-4.5e-5 rad)
        assert!(nut.dpsi.abs() > 1e-8 && nut.dpsi.abs() < 1e-3);
        assert!(nut.deps.abs() > 1e-8 && nut.deps.abs() < 1e-3);
    }
}
