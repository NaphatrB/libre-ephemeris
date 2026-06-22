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

/// Gravitational light deflection by the Sun.
///
/// Based on the relativistic deflection formula:
/// Δ = (1 + γ) * GM / c² * (R_sun / |R_sun|²) × ((R_body - R_sun) / |R_body - R_sun|)
///
/// Where γ = 1 in General Relativity.
/// Implementation follows Meeus Ch. 23 and IAU Resolution B1.3 (2000).
use crate::constants;

/// Sun's gravitational constant in AU³/day²
const GM_SUN: f64 = 2.959122082855911e-4;

/// Compute the gravitational light deflection correction.
///
/// Arguments:
/// - `body_pos`: barycentric position of the observed body (AU)
/// - `earth_pos`: barycentric position of the observer (AU)
/// - `sun_pos`: barycentric position of the Sun (AU) — typically [0;3] for barycentric
///
/// Returns the deflection correction vector to add to the apparent position.
pub fn light_deflection(
    body_pos: &[f64; 3],
    earth_pos: &[f64; 3],
    sun_pos: &[f64; 3],
) -> [f64; 3] {
    let c2 = constants::LE_CLIGHT * constants::LE_CLIGHT;

    // Vector from Sun to Earth (observer)
    let r = [
        earth_pos[0] - sun_pos[0],
        earth_pos[1] - sun_pos[1],
        earth_pos[2] - sun_pos[2],
    ];
    let r_mag = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();

    // Vector from Sun to body
    let r_body = [
        body_pos[0] - sun_pos[0],
        body_pos[1] - sun_pos[1],
        body_pos[2] - sun_pos[2],
    ];
    let r_body_mag = (r_body[0] * r_body[0] + r_body[1] * r_body[1] + r_body[2] * r_body[2]).sqrt();

    // Sun's gravity scale factor
    let g = 2.0 * GM_SUN / (c2 * r_mag);

    // Direction to body from Sun (unit vector)
    let e_body = [
        r_body[0] / r_body_mag,
        r_body[1] / r_body_mag,
        r_body[2] / r_body_mag,
    ];

    // Direction to Earth from Sun (unit vector)
    let e_earth = [r[0] / r_mag, r[1] / r_mag, r[2] / r_mag];

    // Cosine of angle between body and Sun directions
    let cos_theta = e_body[0] * e_earth[0] + e_body[1] * e_earth[1] + e_body[2] * e_earth[2];
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-30);

    // Only apply if the body is not too close to the Sun (avoid singularity)
    if sin_theta < 1e-10 {
        return [0.0; 3];
    }

    // tan(theta/2) factor
    let _tan_half = (1.0 - cos_theta) / sin_theta;

    // Unit vector perpendicular to Sun-Earth line, in the plane of Sun-Earth-body
    // e_⊥ = (e_body - cos_theta * e_earth) / sin_theta
    let e_perp = [
        (e_body[0] - cos_theta * e_earth[0]) / sin_theta,
        (e_body[1] - cos_theta * e_earth[1]) / sin_theta,
        (e_body[2] - cos_theta * e_earth[2]) / sin_theta,
    ];

    let magnitude = g * (1.0 + cos_theta) / sin_theta;

    [
        magnitude * e_perp[0],
        magnitude * e_perp[1],
        magnitude * e_perp[2],
    ]
}
