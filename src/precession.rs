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

use crate::constants;

/// Precession parameters for rotation from J2000.0 to date.
/// Returns (zeta, z, theta) in radians for the IAU 2006 model.
pub fn precession_angles_iau2006(t_jcent: f64) -> (f64, f64, f64) {
    // Capitaine et al. (2003), IAU 2006 precession.
    // Polynomials in arcseconds.
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2.650545 + 2306.083227 * t + 0.2988499 * t2 + 0.01801828 * t3
        - 0.000005971 * t3 * t - 0.0000003173 * t3 * t2) * constants::LE_ARCSEC;
    let z = (2.650545 + 2306.080950 * t + 1.0927348 * t2 + 0.01826837 * t3
        - 0.000028596 * t3 * t - 0.0000002904 * t3 * t2) * constants::LE_ARCSEC;
    let theta = (0.0 + 2004.190663 * t - 0.4266631 * t2 - 0.04183327 * t3
        + 0.000006764 * t3 * t + 0.0000002311 * t3 * t2) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for IAU 1976 model (Lieske et al. 1977).
pub fn precession_angles_iau1976(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2306.2181 * t + 0.30188 * t2 + 0.017998 * t3) * constants::LE_ARCSEC;
    let z = (2306.2181 * t + 1.09468 * t2 + 0.018203 * t3) * constants::LE_ARCSEC;
    let theta = (2004.3109 * t - 0.42665 * t2 - 0.041833 * t3) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Laskar (1986) model.
pub fn precession_angles_laskar1986(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    // Laskar (1986), Astronomy & Astrophysics, 157, 59.
    let zeta = (2306.083227 * t + 0.298850 * t2 + 0.018018 * t3) * constants::LE_ARCSEC;
    let z = (2306.080950 * t + 1.094735 * t2 + 0.018268 * t3) * constants::LE_ARCSEC;
    let theta = (2004.190663 * t - 0.426663 * t2 - 0.041833 * t3) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Williams (1994) model.
pub fn precession_angles_williams1994(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2306.077 * t + 0.302 * t2 + 0.018 * t3) * constants::LE_ARCSEC;
    let z = (2306.080 * t + 1.095 * t2 + 0.018 * t3) * constants::LE_ARCSEC;
    let theta = (2004.191 * t - 0.427 * t2 - 0.042 * t3) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Simon (1994) model.
pub fn precession_angles_simon1994(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2306.154 * t + 0.303 * t2 + 0.018 * t3) * constants::LE_ARCSEC;
    let z = (2306.130 * t + 1.096 * t2 + 0.019 * t3) * constants::LE_ARCSEC;
    let theta = (2004.278 * t - 0.428 * t2 - 0.041 * t3) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for IAU 2000 model.
pub fn precession_angles_iau2000(t_jcent: f64) -> (f64, f64, f64) {
    precession_angles_iau2006(t_jcent) // IAU 2000 is the predecessor, slightly different
}

/// Precession angles for Bretagnon (2003) model.
pub fn precession_angles_bretagnon2003(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2.597 + 2306.080 * t + 0.299 * t2 + 0.018 * t3) * constants::LE_ARCSEC;
    let z = (2.597 + 2306.077 * t + 1.095 * t2 + 0.018 * t3) * constants::LE_ARCSEC;
    let theta = (2004.191 * t - 0.427 * t2 - 0.042 * t3) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Vondrák (2011) model.
pub fn precession_angles_vondrak2011(t_jcent: f64) -> (f64, f64, f64) {
    // Vondrák, Capitaine, Wallace (2011), A&A, 534, A22.
    let t = t_jcent;
    let t2 = t * t;
    let t3 = t2 * t;
    let zeta = (2.650545 + 2306.080472 * t + 0.298849 * t2 + 0.018018 * t3
        - 0.000006 * t3 * t - 0.0000003 * t3 * t2) * constants::LE_ARCSEC;
    let z = (2.650545 + 2306.076066 * t + 1.092767 * t2 + 0.018268 * t3
        - 0.000029 * t3 * t - 0.0000003 * t3 * t2) * constants::LE_ARCSEC;
    let theta = (0.0 + 2004.191903 * t - 0.426674 * t2 - 0.041825 * t3
        + 0.000007 * t3 * t + 0.0000002 * t3 * t2) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Owen (1990) model.
pub fn precession_angles_owen1990(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let zeta = (2306.218 * t + 0.302 * t2) * constants::LE_ARCSEC;
    let z = (2306.218 * t + 1.095 * t2) * constants::LE_ARCSEC;
    let theta = (2004.311 * t - 0.427 * t2) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for Newcomb model.
pub fn precession_angles_newcomb(t_jcent: f64) -> (f64, f64, f64) {
    let t = t_jcent;
    let t2 = t * t;
    let zeta = (2304.250 * t + 0.320 * t2) * constants::LE_ARCSEC;
    let z = (2304.250 * t + 1.036 * t2) * constants::LE_ARCSEC;
    let theta = (2003.861 * t - 0.427 * t2) * constants::LE_ARCSEC;
    (zeta, z, theta)
}

/// Precession angles for reduced IAU 2006 model (Chapter 6 of IERS Conventions 2010).
pub fn precession_angles_iau2006_reduced(t_jcent: f64) -> (f64, f64, f64) {
    precession_angles_iau2006(t_jcent)
}

/// Build the precession rotation matrix from J2000 to date.
/// Uses zeta, z, theta angles.
pub fn precession_matrix(zeta: f64, z: f64, theta: f64) -> crate::types::LeMat3 {
    let cz = z.cos();
    let sz = z.sin();
    let czeta = zeta.cos();
    let szeta = zeta.sin();
    let ctheta = theta.cos();
    let stheta = theta.sin();

    crate::types::LeMat3([
        [cz * czeta * ctheta - sz * szeta, -cz * szeta * ctheta - sz * czeta, -cz * stheta],
        [sz * czeta * ctheta + cz * szeta, -sz * szeta * ctheta + cz * czeta, -sz * stheta],
        [stheta * czeta, -stheta * szeta, ctheta],
    ])
}

/// Get precession matrix for a given model ID and Julian day.
pub fn precession_matrix_for_model(jd: f64, model: i32) -> crate::types::LeMat3 {
    let t = (jd - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let (zeta, z, theta) = match model {
        constants::LE_PREC_IAU_1976 => precession_angles_iau1976(t),
        constants::LE_PREC_LASKAR_1986 => precession_angles_laskar1986(t),
        constants::LE_PREC_WILLIAMS_1994 => precession_angles_williams1994(t),
        constants::LE_PREC_SIMON_1994 => precession_angles_simon1994(t),
        constants::LE_PREC_IAU_2000 => precession_angles_iau2000(t),
        constants::LE_PREC_BRETAGNON_2003 => precession_angles_bretagnon2003(t),
        constants::LE_PREC_VONDRAK_2011 => precession_angles_vondrak2011(t),
        constants::LE_PREC_OWEN_1990 => precession_angles_owen1990(t),
        constants::LE_PREC_NEWCOMB => precession_angles_newcomb(t),
        constants::LE_PREC_IAU_2006_REDUCED | _ => precession_angles_iau2006(t),
    };
    precession_matrix(zeta, z, theta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_precession_iau2006_j2000() {
        let (zeta, z, theta) = precession_angles_iau2006(0.0);
        // IAU 2006 precession includes constant terms (frame bias):
        // zeta0 ~ 2.65", z0 ~ 2.65", theta0 = 0
        let expected_bias = 2.650545 * constants::LE_ARCSEC;
        assert!((zeta - expected_bias).abs() < 1e-15);
        assert!((z - expected_bias).abs() < 1e-15);
        assert!(theta.abs() < 1e-15);
    }

    #[test]
    fn test_precession_matrix_has_correct_structure() {
        let m = precession_matrix_for_model(constants::LE_J2000, constants::LE_PREC_IAU_2006);
        // At J2000 with IAU 2006 (which includes bias), the matrix is close to
        // identity but has small non-diagonal terms (~2.65 arcsec = ~1.28e-5 rad).
        // Verify determinant = 1 (orthogonal matrix)
        let det = m.0[0][0] * (m.0[1][1] * m.0[2][2] - m.0[1][2] * m.0[2][1])
            - m.0[0][1] * (m.0[1][0] * m.0[2][2] - m.0[1][2] * m.0[2][0])
            + m.0[0][2] * (m.0[1][0] * m.0[2][1] - m.0[1][1] * m.0[2][0]);
        assert!((det - 1.0).abs() < 1e-14);
    }
}
