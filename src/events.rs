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
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

/// Event search: find aspects, ingresses, and transits between celestial bodies.
///
/// Uses binary search and Newton's method to find the exact time of
/// planetary events within a given date range.

use crate::calc;
use crate::constants;

/// Aspect type (angular separation between two bodies).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Aspect {
    Conjunction,   // 0°
    Opposition,    // 180°
    Trine,         // 120°
    Square,        // 90°
    Sextile,       // 60°
    SemiSextile,   // 30°
    Quincunx,      // 150°
    SemiSquare,    // 45°
    Sesquiquadrate, // 135°
    Quintile,      // 72°
    BiQuintile,    // 144°
}

impl Aspect {
    pub fn angle_deg(&self) -> f64 {
        match self {
            Aspect::Conjunction => 0.0,
            Aspect::Opposition => 180.0,
            Aspect::Trine => 120.0,
            Aspect::Square => 90.0,
            Aspect::Sextile => 60.0,
            Aspect::SemiSextile => 30.0,
            Aspect::Quincunx => 150.0,
            Aspect::SemiSquare => 45.0,
            Aspect::Sesquiquadrate => 135.0,
            Aspect::Quintile => 72.0,
            Aspect::BiQuintile => 144.0,
        }
    }

    pub fn orb_deg(&self) -> f64 {
        match self {
            Aspect::Conjunction => 8.0,
            Aspect::Opposition => 8.0,
            Aspect::Trine => 7.0,
            Aspect::Square => 6.0,
            Aspect::Sextile => 5.0,
            Aspect::SemiSextile => 3.0,
            Aspect::Quincunx => 3.0,
            Aspect::SemiSquare => 2.0,
            Aspect::Sesquiquadrate => 2.0,
            Aspect::Quintile => 2.0,
            Aspect::BiQuintile => 2.0,
        }
    }
}

/// An astrological event.
#[derive(Debug, Clone)]
pub struct Event {
    pub jd: f64,
    pub body1: i32,
    pub body2: i32,
    pub aspect: Aspect,
    pub separation_deg: f64,
}

/// Compute the angular separation between two bodies at a given JD.
fn separation_at_jd(jd: f64, body1: i32, body2: i32) -> f64 {
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut xx1 = [0.0_f64; 24];
    let mut xx2 = [0.0_f64; 24];
    let mut serr = [0_i8; 256];
    let rc1 = unsafe { calc::le_calc_ut(jd, body1, flags, xx1.as_mut_ptr(), serr.as_mut_ptr()) };
    let rc2 = unsafe { calc::le_calc_ut(jd, body2, flags, xx2.as_mut_ptr(), serr.as_mut_ptr()) };
    if rc1 != 0 || rc2 != 0 {
        return 999.0;
    }
    let dot = xx1[0] * xx2[0] + xx1[1] * xx2[1] + xx1[2] * xx2[2];
    let d1 = (xx1[0] * xx1[0] + xx1[1] * xx1[1] + xx1[2] * xx1[2]).sqrt();
    let d2 = (xx2[0] * xx2[0] + xx2[1] * xx2[1] + xx2[2] * xx2[2]).sqrt();
    if d1 < 1e-15 || d2 < 1e-15 {
        return 999.0;
    }
    (dot / (d1 * d2)).clamp(-1.0, 1.0).acos() * constants::LE_RAD
}

/// Find the exact time of an aspect between two bodies using binary search.
///
/// Arguments:
/// - `jd_start`: start of search range (JD)
/// - `jd_end`: end of search range (JD)
/// - `body1`: first body index
/// - `body2`: second body index
/// - `aspect`: aspect to search for
/// - `tolerance_days`: precision of result (default 1e-4 ≈ 8.6 seconds)
///
/// Returns: the Event if found, or None.
pub fn find_aspect(
    jd_start: f64,
    jd_end: f64,
    body1: i32,
    body2: i32,
    aspect: Aspect,
    tolerance_days: f64,
) -> Option<Event> {
    let target_deg = aspect.angle_deg();
    let target_rad = target_deg * constants::LE_DEG;
    let tol = tolerance_days;

    let sep_start = separation_at_jd(jd_start, body1, body2);
    let sep_end = separation_at_jd(jd_end, body1, body2);

    let diff_start = sep_start - target_rad;
    let diff_end = sep_end - target_rad;

    if diff_start * diff_end > 0.0 {
        return None;
    }

    let mut lo = jd_start;
    let mut hi = jd_end;
    let mut mid = 0.0;

    for _ in 0..60 {
        mid = (lo + hi) / 2.0;
        if hi - lo < tol {
            break;
        }
        let sep_mid = separation_at_jd(mid, body1, body2);
        let diff_mid = sep_mid - target_rad;
        if diff_start * diff_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some(Event {
        jd: mid,
        body1,
        body2,
        aspect,
        separation_deg: separation_at_jd(mid, body1, body2) * constants::LE_RAD,
    })
}

/// Find all aspects between two bodies within a date range.
pub fn find_all_aspects(
    jd_start: f64,
    jd_end: f64,
    body1: i32,
    body2: i32,
    aspects: &[Aspect],
    tolerance_days: f64,
) -> Vec<Event> {
    let mut events = Vec::new();
    for aspect in aspects {
        if let Some(event) = find_aspect(jd_start, jd_end, body1, body2, *aspect, tolerance_days) {
            events.push(event);
        }
    }
    events
}

/// Find the time of an ingress (body entering a sign).
///
/// Returns the JD when the body's ecliptic longitude crosses
/// the given boundary (0°, 30°, 60°, etc.).
pub fn find_ingress(
    jd_start: f64,
    jd_end: f64,
    body: i32,
    sign_boundary_deg: f64,
    tolerance_days: f64,
) -> Option<f64> {
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_ECLIPTIC | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;

    let lon_at = |jd: f64| -> f64 {
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { calc::le_calc_ut(jd, body, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { return 999.0; }
        let lon = xx[1].atan2(xx[0]);
        let mut lon_deg = lon * constants::LE_RAD;
        if lon_deg < 0.0 { lon_deg += 360.0; }
        lon_deg
    };

    let lon_start = lon_at(jd_start);
    let lon_end = lon_at(jd_end);

    // Check if the boundary was crossed: the longitude should pass through
    // the boundary value. Use a simple sign check on the difference.
    let d_start = (lon_start - sign_boundary_deg + 360.0) % 360.0;
    let d_end = (lon_end - sign_boundary_deg + 360.0) % 360.0;
    if d_start * d_end > 0.0 && (d_end - d_start).abs() < 180.0 {
        return None;
    }

    let mut lo = jd_start;
    let mut hi = jd_end;
    let mut mid = 0.0;

    for _ in 0..60 {
        mid = (lo + hi) / 2.0;
        if hi - lo < tolerance_days {
            break;
        }
        let lon_mid = lon_at(mid);
        let diff_mid = (lon_mid - sign_boundary_deg + 360.0) % 360.0;
        if d_start * diff_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some(mid)
}

/// Find the time of a transit (body crossing the meridian).
///
/// Returns the JD when the body's RA equals the local sidereal time.
pub fn find_transit(
    jd_start: f64,
    jd_end: f64,
    body: i32,
    longitude: f64,
    tolerance_days: f64,
) -> Option<f64> {
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;

    let ha_at = |jd: f64| -> f64 {
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { calc::le_calc_ut(jd, body, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { return 999.0; }
        let ra = xx[1].atan2(xx[0]);
        let gmst = crate::riseset::mean_sidereal_time_greenwich(jd);
        let lst = gmst + longitude;
        let ha = lst - ra;
        ha.rem_euclid(2.0 * std::f64::consts::PI)
    };

    let ha_start = ha_at(jd_start);
    let ha_end = ha_at(jd_end);

    if ha_start * ha_end > 0.0 {
        return None;
    }

    let mut lo = jd_start;
    let mut hi = jd_end;
    let mut mid = 0.0;

    for _ in 0..60 {
        mid = (lo + hi) / 2.0;
        if hi - lo < tolerance_days {
            break;
        }
        let ha_mid = ha_at(mid);
        if ha_start * ha_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some(mid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aspect_angles() {
        assert!((Aspect::Conjunction.angle_deg() - 0.0).abs() < 1e-10);
        assert!((Aspect::Opposition.angle_deg() - 180.0).abs() < 1e-10);
        assert!((Aspect::Trine.angle_deg() - 120.0).abs() < 1e-10);
    }

    #[test]
    fn test_separation_at_jd() {
        let sep = separation_at_jd(2451545.0, constants::LE_SUN, constants::LE_MOON);
        assert!(sep >= 0.0 && sep < 180.0, "separation={}", sep);
    }
}
