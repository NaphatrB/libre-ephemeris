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

/// Atmospheric refraction and horizontal coordinate calculations.
///
/// Computes apparent altitude from true altitude (refraction),
/// and converts between equatorial (RA/Dec) and horizontal (alt/az) coordinates.

use crate::constants;

/// Convert equatorial coordinates (RA, Dec) to horizontal coordinates (alt, az).
///
/// Azimuth measured from south through west (0° = south, 90° = west).
/// Meeus (1998) Ch. 13.
pub fn equatorial_to_horizontal(ra: f64, dec: f64, lst: f64, lat: f64) -> (f64, f64) {
    let ha = lst - ra;
    let (s_ha, c_ha) = ha.sin_cos();
    let (s_lat, c_lat) = lat.sin_cos();
    let (s_dec, c_dec) = dec.sin_cos();

    let alt = (s_lat * s_dec + c_lat * c_dec * c_ha).asin();
    let az = s_ha.atan2(c_ha * s_lat - s_dec / c_dec * c_lat);
    let az = if az < 0.0 { az + 2.0 * std::f64::consts::PI } else { az };
    (alt, az)
}

/// Convert horizontal coordinates (alt, az) to equatorial (RA, Dec).
///
/// Azimuth measured from south through west (0° = south, 90° = west).
/// Meeus (1998) Ch. 13.
pub fn horizontal_to_equatorial(alt: f64, az: f64, lst: f64, lat: f64) -> (f64, f64) {
    let (s_alt, c_alt) = alt.sin_cos();
    let (s_az, c_az) = az.sin_cos();
    let (s_lat, c_lat) = lat.sin_cos();

    let dec = (s_lat * s_alt - c_lat * c_alt * c_az).asin();
    let ha = s_az.atan2(c_az * s_lat + s_alt / c_alt * c_lat);
    let ra = lst - ha;
    (ra, dec)
}

/// Compute atmospheric refraction correction.
///
/// Bennett (1982) formula: R = cot(h + 7.31/(h + 4.4))
/// where h is altitude in degrees, R is in arcminutes.
/// Returns the correction in radians.
///
/// Arguments:
/// - `alt_true`: true altitude in radians
/// - `pressure`: atmospheric pressure in millibars (default 1010)
/// - `temperature`: temperature in Celsius (default 10)
pub fn refraction(alt_true: f64, pressure: f64, temperature: f64) -> f64 {
    if alt_true <= -0.1 {
        return 0.0;
    }
    let h_deg = alt_true * constants::LE_RAD;
    let arg_deg = h_deg + 7.31 / (h_deg + 4.4);
    let r_arcmin = 1.0 / arg_deg.to_radians().tan();
    let r_arcmin = r_arcmin.clamp(0.0, 60.0);
    let p_factor = pressure / 1010.0;
    let t_factor = 283.0 / (273.0 + temperature);
    r_arcmin * p_factor * t_factor * constants::LE_DEG / 60.0
}

/// Compute apparent altitude from true altitude (includes refraction).
pub fn apparent_altitude(alt_true: f64, pressure: f64, temperature: f64) -> f64 {
    alt_true + refraction(alt_true, pressure, temperature)
}

/// Compute true altitude from apparent altitude (removes refraction).
pub fn true_altitude(alt_apparent: f64, pressure: f64, temperature: f64) -> f64 {
    let mut alt = alt_apparent;
    for _ in 0..3 {
        let r = refraction(alt, pressure, temperature);
        alt = alt_apparent - r;
    }
    alt
}

/// Compute azimuth of a body from its equatorial coordinates.
pub fn azimuth(ra: f64, dec: f64, lst: f64, lat: f64) -> f64 {
    let (_alt, az) = equatorial_to_horizontal(ra, dec, lst, lat);
    az
}

/// Compute altitude of a body from its equatorial coordinates.
pub fn altitude(ra: f64, dec: f64, lst: f64, lat: f64) -> f64 {
    let (alt, _az) = equatorial_to_horizontal(ra, dec, lst, lat);
    alt
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_refraction_positive() {
        let r = refraction(30.0 * constants::LE_DEG, 1010.0, 10.0);
        assert!(r > 0.0 && r < 0.1, "refraction={:.6} rad", r);
    }

    #[test]
    fn test_refraction_zero_at_horizon() {
        let r = refraction(0.0, 1010.0, 10.0);
        assert!(r > 0.0, "refraction at horizon should be positive");
    }

    #[test]
    fn test_equatorial_to_horizontal_roundtrip() {
        let ra = 5.0;
        let dec = 0.5;
        let lst = 6.0;
        let lat = 0.7;
        let (alt, az) = equatorial_to_horizontal(ra, dec, lst, lat);
        let (ra2, dec2) = horizontal_to_equatorial(alt, az, lst, lat);
        let ra_diff = (ra - ra2 + 2.0 * std::f64::consts::PI) % (2.0 * std::f64::consts::PI);
        assert!(ra_diff < 1e-10 || (2.0 * std::f64::consts::PI - ra_diff) < 1e-10,
            "ra diff={}", ra_diff);
        assert!((dec - dec2).abs() < 1e-10, "dec diff={}", dec - dec2);
    }
}
