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
use crate::types::LeMat3;

/// Frame bias matrix from ICRS to FK5 J2000 (or vice versa).
///
/// The ICRS (International Celestial Reference System) is aligned with the
/// extragalactic radio frame. The FK5 (J2000.0) is the previous optical
/// reference frame. The difference is a few milliarcseconds.
///
/// IAU 2006 Resolution B2 defines the bias parameters:
/// ξ0 = -16.6170 mas  (x-axis offset)
/// η0 = -6.8192 mas   (y-axis offset)
/// da0 = -14.6 mas    (equator offset)
///
/// Source: IERS Conventions 2010, Chapter 5.
const BIAS_XI_0: f64 = -16.6170;    // mas
const BIAS_ETA_0: f64 = -6.8192;     // mas
const BIAS_DA_0: f64 = -14.6;        // mas

fn mas_to_rad(mas: f64) -> f64 {
    mas * 4.848136811095359e-9 // 1 mas = π/(180*3600*1000)
}

/// Build the frame bias matrix (ICRS -> FK5 J2000).
/// This is the standard IAU 2006 bias matrix using small-angle approximation.
///
/// The three bias parameters are:
///   ξ0 = −16.617 mas (rotation about y-axis, toward the ICRS pole)
///   η0 = −6.8192 mas  (rotation about x-axis, along the equinox)
///   da0 = −14.6 mas   (rotation about z-axis, equator offset)
///
/// Bias = R_z(-da0) * R_y(ξ0) * R_x(η0)
/// For milliarcsecond angles, the small-angle approximation is exact.
pub fn frame_bias_matrix() -> LeMat3 {
    let xi0 = mas_to_rad(BIAS_XI_0);
    let eta0 = mas_to_rad(BIAS_ETA_0);
    let da0 = mas_to_rad(BIAS_DA_0);

    // Small-angle rotation matrix (valid for angles << 1 rad):
    //   R_z(-da) * R_y(xi) * R_x(eta) ≈
    //   [ 1,   -da,   xi ]
    //   [ da,   1,   -eta ]
    //   [-xi,   eta,   1  ]
    LeMat3([
        [1.0, -da0, xi0],
        [da0, 1.0, -eta0],
        [-xi0, eta0, 1.0],
    ])
}

/// Apply frame bias to a position vector (ICRS -> FK5).
pub fn apply_frame_bias(pos: &[f64; 3], bias_model: i32) -> [f64; 3] {
    match bias_model {
        constants::LE_BIAS_NONE => *pos,
        constants::LE_BIAS_IAU_2000 | constants::LE_BIAS_IAU_2006 => {
            frame_bias_matrix().transform(pos)
        }
        _ => *pos,
    }
}
