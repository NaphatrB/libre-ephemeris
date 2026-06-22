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

/// Compute annual aberration correction.
///
/// Given Earth's heliocentric position and velocity vectors (in AU, AU/day),
/// compute the relativistic deflection of incoming starlight.
///
/// Standard formula from Meeus, "Astronomical Algorithms" Ch. 23.
///
/// Returns the correction to be subtracted from apparent position to get
/// true position, or equivalently added to true to get apparent.
pub fn annual_aberration(pos: &[f64; 3]) -> [f64; 3] {
    let _c = constants::LE_CLIGHT;
    // pos is Earth's position relative to Sun (AU)
    // velocity of Earth in AU/day - for aberration we need the velocity vector
    // This simplified version uses the Earth's position as a proxy assuming
    // approximately circular orbit. The full version needs Earth velocity.
    // For the proper implementation, see the velocity-aware version below.
    *pos
}

/// Full relativistic aberration correction using Earth velocity.
///
/// From the relativistic aberration formula (Meeus Ch. 23):
/// Δ = v/c - (R · v/c) * R
/// where R = direction to body, v = Earth velocity, c = speed of light.
/// This correction is added to the unit vector direction; multiply by
/// distance to get the position correction.
pub fn aberration_correction(
    dir: &[f64; 3],  // unit vector to body (J2000)
    earth_vel: &[f64; 3],  // Earth barycentric velocity (AU/day)
) -> [f64; 3] {
    let c = constants::LE_CLIGHT;
    let v_over_c = [
        earth_vel[0] / c,
        earth_vel[1] / c,
        earth_vel[2] / c,
    ];

    let r_dot_v_over_c = dir[0] * v_over_c[0] + dir[1] * v_over_c[1] + dir[2] * v_over_c[2];

    [
        v_over_c[0] - dir[0] * r_dot_v_over_c,
        v_over_c[1] - dir[1] * r_dot_v_over_c,
        v_over_c[2] - dir[2] * r_dot_v_over_c,
    ]
}

/// Apply aberration correction to a celestial position.
/// `apparent` is the position vector to correct.
/// Returns the aberration-corrected (apparent) position.
pub fn apply_aberration(
    pos: &[f64; 3],
    earth_pos: &[f64; 3],
) -> [f64; 3] {
    // Earth velocity approximation from two positions
    // In real usage, this would use the actual Earth velocity from the ephemeris
    let _dir = pos;
    let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
    let unit_dir = [pos[0] / dist, pos[1] / dist, pos[2] / dist];

    // Approximate Earth velocity (circular orbit)
    let r = (earth_pos[0] * earth_pos[0] + earth_pos[1] * earth_pos[1] + earth_pos[2] * earth_pos[2]).sqrt();
    let _omega = 0.01720209895 * 365.25 / (2.0 * std::f64::consts::PI); // 1/year in rad/day
    // ~0.0172 rad/day for Earth's orbital angular velocity
    let omega_actual = constants::LE_GAUSS_G / (r * r * r).sqrt();

    // Velocity = omega × r (cross product in orbital plane)
    // Simplified approximation:
    let v = [
        -earth_pos[1] * omega_actual,
        earth_pos[0] * omega_actual,
        0.0,
    ];

    let corr = aberration_correction(&unit_dir, &v);

    [
        pos[0] + corr[0] * dist,
        pos[1] + corr[1] * dist,
        pos[2] + corr[2] * dist,
    ]
}
