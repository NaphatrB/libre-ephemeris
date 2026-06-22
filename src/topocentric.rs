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

/// Earth's equatorial radius in km (WGS84)
const EARTH_RADIUS: f64 = 6378.137;

/// Earth's flattening (WGS84)
const EARTH_FLATTENING: f64 = 1.0 / 298.257223563;

/// Convert geodetic coordinates to geocentric equatorial position vector (in AU).
///
/// Uses WGS84 ellipsoid. Returns observer's position in Earth-centered inertial frame,
/// which can be combined with Earth's barycentric position for topocentric corrections.
///
/// Arguments:
/// - `lon`: geographic longitude in degrees (positive east)
/// - `lat`: geographic latitude in degrees (positive north)
/// - `alt`: altitude above WGS84 ellipsoid in meters
/// - `lst`: local sidereal time in radians
///
/// Returns: observer position vector in AU (Earth-centered, equatorial J2000)
pub fn geodetic_to_geocentric(
    lon: f64,
    lat: f64,
    alt: f64,
    lst: f64,
) -> [f64; 3] {
    let phi = lat * constants::LE_DEG;
    let _lambda = lon * constants::LE_DEG;

    let sin_phi = phi.sin();
    let cos_phi = phi.cos();

    let ecc2 = 2.0 * EARTH_FLATTENING - EARTH_FLATTENING * EARTH_FLATTENING;
    let denom = (1.0 - ecc2 * sin_phi * sin_phi).sqrt();

    // Distance from Earth's center to observer (in km)
    let c = EARTH_RADIUS / denom + alt / 1000.0;
    let s = EARTH_RADIUS * (1.0 - ecc2) / denom + alt / 1000.0;

    let x_km = c * cos_phi * lst.cos();
    let y_km = c * cos_phi * lst.sin();
    let z_km = s * sin_phi;

    // Convert km to AU
    let km_to_au = 1.0 / 1.495978707e8;
    [x_km * km_to_au, y_km * km_to_au, z_km * km_to_au]
}

/// Compute topocentric correction for a body position.
///
/// Given the observer's geocentric position (from `geodetic_to_geocentric`)
/// and the Earth's barycentric position (from ephemeris), compute the
/// topocentric position of a body.
///
/// The topocentric position is:
/// pos_topocentric = pos_geocentric - obs_geocentric
/// where all vectors are relative to Earth's center.
///
/// For a body at distance d, the correction is primarily the parallax:
/// π = arcsin(R_earth / d) where d is the body's distance.
pub fn topocentric_correction(
    body_earth_pos: &[f64; 3],  // body position relative to Earth (AU)
    obs_geocentric: &[f64; 3],  // observer position relative to Earth center (AU)
) -> [f64; 3] {
    [
        body_earth_pos[0] - obs_geocentric[0],
        body_earth_pos[1] - obs_geocentric[1],
        body_earth_pos[2] - obs_geocentric[2],
    ]
}

/// Compute diurnal parallax (horizontal parallax) in radians.
///
/// For a body at distance d (AU) from Earth's center,
/// the equatorial horizontal parallax is:
/// π = arcsin(R_earth / d)
pub fn horizontal_parallax(distance_au: f64) -> f64 {
    let r_earth_au = EARTH_RADIUS / 1.495978707e8;
    (r_earth_au / distance_au).asin()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geodetic_at_equator() {
        // Observer at equator, lon=0, lat=0, sea level
        let pos = geodetic_to_geocentric(0.0, 0.0, 0.0, 0.0);
        // Should be on x-axis
        assert!(pos[1].abs() < 1e-15);
        assert!(pos[2].abs() < 1e-15);
        assert!(pos[0] > 0.0);
    }

    #[test]
    fn test_parallax_moon() {
        // Moon distance ~0.0025 AU, parallax should be ~1 degree
        let p = horizontal_parallax(0.0025);
        let p_deg = p / constants::LE_DEG;
        assert!((p_deg - 0.95).abs() < 0.1);
    }
}
