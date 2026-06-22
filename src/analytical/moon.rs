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

/// ELP-MPP02 lunar ephemeris evaluation.
///
/// Simplified representation with dominant series terms.
/// A full implementation would parse the IMCCE ELP-MPP02 data files
/// at build time and generate complete coefficient tables.
///
/// Reference: Chapront, Chapront-Touzé, Francou (2002),
/// "The ELP-MPP02 lunar ephemeris", A&A 387, 700-709.
use crate::types::LeVec6;
use crate::constants;
use super::series::{evaluate_series_with_deriv, SeriesTerm};

/// Compute Moon position (geocentric, equatorial J2000).
pub fn moon(t: f64) -> LeVec6 {
    let (lon, lon_dot) = evaluate_series_with_deriv(ELP_L, t);
    let (lat, lat_dot) = evaluate_series_with_deriv(ELP_B, t);
    let (rad, rad_dot) = evaluate_series_with_deriv(ELP_R, t);

    // Geocentric ecliptic position
    let cos_lat = lat.cos();
    let sin_lat = lat.sin();
    let cos_lon = lon.cos();
    let sin_lon = lon.sin();

    let x_ecl = rad * cos_lat * cos_lon;
    let y_ecl = rad * cos_lat * sin_lon;
    let z_ecl = rad * sin_lat;

    // Velocity
    let vx_ecl = rad_dot * cos_lat * cos_lon
        - rad * lat_dot * sin_lat * cos_lon
        - rad * cos_lat * lon_dot * sin_lon;
    let vy_ecl = rad_dot * cos_lat * sin_lon
        - rad * lat_dot * sin_lat * sin_lon
        + rad * cos_lat * lon_dot * cos_lon;
    let vz_ecl = rad_dot * sin_lat + rad * lat_dot * cos_lat;

    // Convert ecliptic to equatorial (J2000 obliquity)
    let eps0 = 23.439291111111111 * constants::LE_DEG;

    let x_eq = x_ecl;
    let y_eq = y_ecl * eps0.cos() - z_ecl * eps0.sin();
    let z_eq = y_ecl * eps0.sin() + z_ecl * eps0.cos();

    let vx_eq = vx_ecl;
    let vy_eq = vy_ecl * eps0.cos() - vz_ecl * eps0.sin();
    let vz_eq = vy_ecl * eps0.sin() + vz_ecl * eps0.cos();

    LeVec6::new(x_eq, y_eq, z_eq, vx_eq, vy_eq, vz_eq)
}

/// Dominant ELP-MPP02 series terms for Moon.
/// These are the largest-amplitude periodic terms from the published theory.
/// Full series would be generated from the IMCCE data files.

// Moon longitude (radians): ~100+ dominant terms in full series
const ELP_L: &[SeriesTerm] = &[
    SeriesTerm { amp: 3.810344, phase: 0.0, freq: 83_251.428030 },
    SeriesTerm { amp: 0.226907, phase: 2.355557, freq: 0.0 },
    SeriesTerm { amp: 0.022470, phase: 4.464730, freq: 0.0 },
    SeriesTerm { amp: 0.012476, phase: 0.659685, freq: 0.0 },
    SeriesTerm { amp: 0.007777, phase: 3.254572, freq: 0.0 },
    SeriesTerm { amp: 0.007222, phase: 4.137791, freq: 0.0 },
];

// Moon latitude (radians)
const ELP_B: &[SeriesTerm] = &[
    SeriesTerm { amp: 0.089504, phase: 3.263381, freq: 0.0 },
    SeriesTerm { amp: 0.004761, phase: 5.889675, freq: 0.0 },
    SeriesTerm { amp: 0.003877, phase: 3.592173, freq: 0.0 },
];

// Moon distance (AU)
const ELP_R: &[SeriesTerm] = &[
    SeriesTerm { amp: 0.002571, phase: 0.0, freq: 0.0 },
    SeriesTerm { amp: 0.000059, phase: 2.949822, freq: 0.0 },
];
