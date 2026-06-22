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

/// Rise, set, and transit time calculations.
///
/// Given a celestial body's apparent position and an observer's location,
/// computes the UT times of rising, setting, and meridian transit.
///
/// Reference: Meeus, "Astronomical Algorithms" (2nd ed.), Chapter 15.
use crate::constants;
use crate::transform::csnorm;

const DEG_TO_RAD: f64 = constants::LE_DEG;

/// Rise/set/transit event codes (matching `rsmi` parameter conventions).
pub const RISE: i32 = 0;
pub const SET: i32 = 1;
pub const TRANSIT_UPPER: i32 = 2;  // upper culmination (transit)
pub const TRANSIT_LOWER: i32 = 3; // lower culmination

/// Standard horizon altitude (degrees) for rise/set.
/// Refraction at horizon ≈ -0.5667°.
const HORIZON_STD: f64 = -0.5667;

/// Solar horizon: includes refraction + apparent semi-diameter ≈ -0.8333°.
const HORIZON_SUN: f64 = -0.8333;

pub fn mean_sidereal_time_greenwich(jd: f64) -> f64 {
    let t = (jd - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let gmst = 280.46061837 + 360.98564736629 * (jd - constants::LE_J2000)
        + 0.000387933 * t * t - t * t * t / 38710000.0;
    csnorm(gmst)
}

/// Compute the hour angle H (in radians) for a given altitude.
fn hour_angle(altitude_deg: f64, lat_deg: f64, dec_rad: f64) -> Option<f64> {
    let alt = altitude_deg * DEG_TO_RAD;
    let lat = lat_deg * DEG_TO_RAD;
    let sin_alt = alt.sin();
    let sin_lat = lat.sin();
    let cos_lat = lat.cos();
    let sin_dec = dec_rad.sin();
    let cos_dec = dec_rad.cos();

    let cos_h = (sin_alt - sin_lat * sin_dec) / (cos_lat * cos_dec);
    if cos_h < -1.0 || cos_h > 1.0 {
        return None; // Body is circumpolar or never rises
    }
    Some(cos_h.acos())
}

/// Compute rise, set, or transit time for a given body.
///
/// Returns the Julian day (UT) of the event, or an error code.
pub fn rise_trans(
    tjd_start: f64,
    ra_rad: f64,
    dec_rad: f64,
    geo_lon_deg: f64,
    geo_lat_deg: f64,
    rsmi: i32,
    is_sun: bool,
) -> Result<f64, i32> {
    let horizon = if is_sun { HORIZON_SUN } else { HORIZON_STD };

    let ha = match hour_angle(horizon, geo_lat_deg, dec_rad) {
        Some(h) => h,
        None => return Err(constants::ERR_OUT_OF_RANGE),
    };

    let ha_deg = ha / DEG_TO_RAD;

    // Local Sidereal Time at event
    let lst_event_deg = match rsmi {
        RISE => (ra_rad / DEG_TO_RAD - ha_deg).rem_euclid(360.0),
        SET => (ra_rad / DEG_TO_RAD + ha_deg).rem_euclid(360.0),
        TRANSIT_UPPER => (ra_rad / DEG_TO_RAD).rem_euclid(360.0),
        TRANSIT_LOWER => (ra_rad / DEG_TO_RAD + 180.0).rem_euclid(360.0),
        _ => return Err(constants::LE_ERR_INVALID_PARAMS),
    };

    // GST at event = LST - longitude
    let gst_event = csnorm(lst_event_deg - geo_lon_deg);

    // Compute GST at 0h UT of the start date
    let jd0 = tjd_start.floor() - 0.5;
    let gst0 = mean_sidereal_time_greenwich(jd0);

    // Fraction of day = (GST_event - GST_0) / 360.98564736629
    let mut day_frac = (gst_event - gst0).rem_euclid(360.0) / 360.98564736629;
    if day_frac < 0.0 {
        day_frac += 1.0;
    }

    let tjd_event = jd0 + day_frac;

    // Iterate once to improve accuracy
    let ra2 = ra_rad;
    let dec2 = dec_rad;
    let ha2 = match hour_angle(horizon, geo_lat_deg, dec2) {
        Some(h) => h,
        None => return Err(constants::ERR_OUT_OF_RANGE),
    };
    let ha2_deg = ha2 / DEG_TO_RAD;

    let lst_event_deg2 = match rsmi {
        RISE => (ra2 / DEG_TO_RAD - ha2_deg).rem_euclid(360.0),
        SET => (ra2 / DEG_TO_RAD + ha2_deg).rem_euclid(360.0),
        TRANSIT_UPPER => (ra2 / DEG_TO_RAD).rem_euclid(360.0),
        TRANSIT_LOWER => (ra2 / DEG_TO_RAD + 180.0).rem_euclid(360.0),
        _ => return Err(constants::LE_ERR_INVALID_PARAMS),
    };

    let gst_event2 = csnorm(lst_event_deg2 - geo_lon_deg);
    let gst0_2 = mean_sidereal_time_greenwich(tjd_event.floor() - 0.5);
    let mut day_frac2 = (gst_event2 - gst0_2).rem_euclid(360.0) / 360.98564736629;
    if day_frac2 < 0.0 { day_frac2 += 1.0; }

    Ok(tjd_event.floor() - 0.5 + day_frac2)
}

/// Get geocentric equatorial position for a given body from the analytical engine.
fn get_body_xyz(tjd: f64, ipl: i32) -> Result<[f64; 3], i32> {
    #[cfg(feature = "analytical")]
    {
        let raw = crate::analytical::compute_position(tjd, ipl).map_err(|e| e)?;
        Ok([raw.0[0], raw.0[1], raw.0[2]])
    }
    #[cfg(not(feature = "analytical"))]
    {
        let _ = tjd; let _ = ipl;
        Err(constants::ERR_ENGINE)
    }
}

/// Compute geocentric RA/Dec for a given body.
fn get_body_ra_dec(tjd: f64, ipl: i32) -> Result<(f64, f64), i32> {
    let (x, y, z) = if ipl == 0 {
        // Sun: geocentric = -Earth (heliocentric Earth position)
        let e = get_body_xyz(tjd, constants::LE_EARTH)?;
        (-e[0], -e[1], -e[2])
    } else if ipl == constants::LE_MOON {
        // Moon: already geocentric from analytical engine
        let m = get_body_xyz(tjd, constants::LE_MOON)?;
        (m[0], m[1], m[2])
    } else {
        // Planet: geocentric = planet - Earth
        let e = get_body_xyz(tjd, constants::LE_EARTH)?;
        let p = get_body_xyz(tjd, ipl)?;
        (p[0] - e[0], p[1] - e[1], p[2] - e[2])
    };

    let dist = (x * x + y * y + z * z).sqrt();
    if dist < 1e-15 {
        return Err(constants::ERR_ENGINE);
    }
    let ra = y.atan2(x);
    let dec = (z / dist).asin();
    Ok((ra, dec))
}

/// C ABI: compute rise, set, or transit time.
///
/// For fast-moving bodies (Moon), iterates up to 3 times,
/// recomputing the body position at each estimated event time.
#[no_mangle]
pub unsafe extern "C" fn le_rise_trans(
    tjd_ut: f64,
    ipl: i32,
    _iflag: i32,
    rsmi: i32,
    geo_lon: f64,
    geo_lat: f64,
    _attr: *mut f64,
    serr: *mut i8,
) -> f64 {
    let is_sun = ipl == 0;
    let is_moon = ipl == constants::LE_MOON;
    let max_iter = if is_moon { 3 } else { 1 };

    let mut tjd_est = tjd_ut;
    for iter in 0..max_iter {
        let dt_seconds = crate::delta_t::le_deltat(tjd_est);
        let tjd_et = tjd_est + dt_seconds / 86400.0;

        let (ra_rad, dec_rad) = match get_body_ra_dec(tjd_et, ipl) {
            Ok(p) => p,
            Err(e) => return e as f64,
        };

        match rise_trans(tjd_est, ra_rad, dec_rad, geo_lon, geo_lat, rsmi, is_sun) {
            Ok(tjd_event) => {
                if iter + 1 >= max_iter {
                    if !serr.is_null() { unsafe { *serr = 0; } }
                    return tjd_event;
                }
                // Iterate using the computed event time as the new estimate
                tjd_est = tjd_event;
                // Clamp to reasonable range around the original time
                if tjd_est < tjd_ut - 2.0 { tjd_est = tjd_ut - 2.0; }
                if tjd_est > tjd_ut + 2.0 { tjd_est = tjd_ut + 2.0; }
            }
            Err(e) => {
                if !serr.is_null() { unsafe { *serr = 0; } }
                return e as f64;
            }
        }
    }
    // Should not reach here
    constants::ERR_ENGINE as f64
}
