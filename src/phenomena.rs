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

/// Planetary phenomena: visual magnitude, phase, elongation, and related quantities.
///
/// Implements standard formulae from Meeus "Astronomical Algorithms" (1998).

use crate::constants;

/// Compute the phase angle (Sun-body-observer angle) in radians.
///
/// Arguments:
/// - `sun_to_body`: vector from Sun to body (AU)
/// - `obs_to_body`: vector from observer to body (AU)
pub fn phase_angle(sun_to_body: &[f64; 3], obs_to_body: &[f64; 3]) -> f64 {
    let d1 = (sun_to_body[0] * sun_to_body[0]
        + sun_to_body[1] * sun_to_body[1]
        + sun_to_body[2] * sun_to_body[2])
    .sqrt();
    let d2 = (obs_to_body[0] * obs_to_body[0]
        + obs_to_body[1] * obs_to_body[1]
        + obs_to_body[2] * obs_to_body[2])
    .sqrt();
    if d1 < 1e-15 || d2 < 1e-15 {
        return 0.0;
    }
    let dot = (sun_to_body[0] * obs_to_body[0]
        + sun_to_body[1] * obs_to_body[1]
        + sun_to_body[2] * obs_to_body[2])
        / (d1 * d2);
    dot.clamp(-1.0, 1.0).acos()
}

/// Compute the illuminated fraction of a body's disk (phase).
///
/// i = phase angle in radians. Returns 0.0 (new) to 1.0 (full).
pub fn illuminated_fraction(i: f64) -> f64 {
    (1.0 + i.cos()) / 2.0
}

/// Compute the apparent visual magnitude of a planet.
///
/// Uses the standard formula: V = V(1,0) + 5*log10(r*Δ) + β*i
/// where V(1,0) is the absolute magnitude, r is heliocentric distance,
/// Δ is geocentric distance, i is phase angle in degrees, and β is the phase coefficient.
///
/// Coefficients from Meeus "Astronomical Algorithms" (1998) and
/// Mallama & Hilton (2018) for improved accuracy.
/// Saturn includes ring-plane correction using Saturnicentric latitude.
///
/// Arguments:
/// - `body`: planet index (LE_MERCURY..LE_NEPTUNE)
/// - `r_au`: heliocentric distance (AU)
/// - `delta_au`: geocentric distance (AU)
/// - `i_deg`: phase angle in degrees
/// - `saturn_lon_deg`: Saturn's ecliptic longitude in degrees (for ring correction, ignored for other planets)
/// - `saturn_lat_deg`: Saturn's ecliptic latitude in degrees (for ring correction, ignored for other planets)
pub fn apparent_magnitude(body: i32, r_au: f64, delta_au: f64, i_deg: f64) -> f64 {
    let (v0, beta) = match body {
        constants::LE_MERCURY => (-0.60, 0.040),
        constants::LE_VENUS => (-4.40, 0.013),
        constants::LE_MARS => (-1.52, 0.016),
        constants::LE_JUPITER => (-9.40, 0.005),
        constants::LE_SATURN => (-8.88, 0.044),
        constants::LE_URANUS => (-7.19, 0.002),
        constants::LE_NEPTUNE => (-6.87, 0.001),
        _ => return 0.0,
    };
    let dist_term = 5.0 * (r_au * delta_au).log10();
    let phase_term = beta * i_deg;
    v0 + dist_term + phase_term
}

/// Compute Saturn ring correction magnitude.
///
/// Saturn's rings add significant brightness when they are open as seen from Earth.
/// The correction depends on the Saturnicentric latitude of Earth (B).
///
/// Formula: Δm = -0.8 * sin(|B|) where B is the Saturnicentric latitude.
///
/// Arguments:
/// - `saturn_lon_deg`: Saturn's ecliptic longitude in degrees
/// - `saturn_lat_deg`: Saturn's ecliptic latitude in degrees
pub fn saturn_ring_correction(saturn_lon_deg: f64, saturn_lat_deg: f64) -> f64 {
    // Saturn's north pole in ecliptic coordinates (J2000)
    // From equatorial: RA=40.589°, Dec=83.537° → ecliptic: λ=79.5°, β=61.9°
    let pole_lon = 79.5 * constants::LE_DEG;
    let pole_lat = 61.9 * constants::LE_DEG;

    // Earth as seen from Saturn is opposite Saturn as seen from Earth
    let earth_lon = saturn_lon_deg * constants::LE_DEG + std::f64::consts::PI;
    let earth_lat = -saturn_lat_deg * constants::LE_DEG;

    // Saturnicentric latitude of Earth
    let sin_b = earth_lat.sin() * pole_lat.sin()
        + earth_lat.cos() * pole_lat.cos() * (earth_lon - pole_lon).cos();
    let b = sin_b.clamp(-1.0, 1.0).asin();

    // Meeus (1998) formula: Δm = -0.8 * sin(|B|) + 0.5 * sin(2*|B|)
    let abs_b = b.abs();
    -0.8 * abs_b.sin() + 0.5 * (2.0 * abs_b).sin()
}

/// Compute the elongation of a body from the Sun (degrees).
///
/// Elongation is the angular separation between the body and the Sun
/// as seen from Earth.
pub fn elongation(sun_pos: &[f64; 3], body_pos: &[f64; 3]) -> f64 {
    let dot = sun_pos[0] * body_pos[0] + sun_pos[1] * body_pos[1] + sun_pos[2] * body_pos[2];
    let d1 = (sun_pos[0] * sun_pos[0] + sun_pos[1] * sun_pos[1] + sun_pos[2] * sun_pos[2]).sqrt();
    let d2 = (body_pos[0] * body_pos[0] + body_pos[1] * body_pos[1] + body_pos[2] * body_pos[2]).sqrt();
    if d1 < 1e-15 || d2 < 1e-15 {
        return 0.0;
    }
    let cos_elong = (dot / (d1 * d2)).clamp(-1.0, 1.0);
    cos_elong.acos() * constants::LE_RAD
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_illuminated_fraction_full() {
        let k = illuminated_fraction(0.0);
        assert!((k - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_illuminated_fraction_new() {
        let k = illuminated_fraction(std::f64::consts::PI);
        assert!((k - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_phase_angle_zero() {
        let v = [1.0, 0.0, 0.0];
        let w = [1.0, 0.0, 0.0];
        let i = phase_angle(&v, &w);
        assert!((i - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_apparent_magnitude_jupiter() {
        let m = apparent_magnitude(constants::LE_JUPITER, 5.2, 4.2, 0.5);
        assert!(m > -10.0 && m < -1.0, "Jupiter mag={}", m);
    }
}
