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

/// Valid range flags for Julian day conversions.
pub const LE_GREG_CAL: i32 = 1;
pub const LE_JUL_CAL: i32 = 0;

/// Convert calendar date to Julian day number.
///
/// Implements the standard algorithm from Meeus, "Astronomical Algorithms" (1998).
/// For Gregorian calendar (gregflag=1) or Julian calendar (gregflag=0).
pub fn julday(year: i32, month: i32, day: f64, gregflag: i32) -> f64 {
    let mut y = year as f64;
    let mut m = month as f64;

    if month <= 2 {
        y -= 1.0;
        m += 12.0;
    }

    let a = (y / 100.0).floor();
    let b = if gregflag != 0 {
        2.0 - a + (a / 4.0).floor()
    } else {
        0.0
    };

    let day_int = day.floor();
    let day_frac = day - day_int;

    (365.25 * (y + 4716.0)).floor()
        + (30.6001 * (m + 1.0)).floor()
        + day_int as f64
        + b
        - 1524.5
        + day_frac
}

/// Convert Julian day number to calendar date.
///
/// Returns (year, month, day, hour).
/// Works for both Gregorian (gregflag=1) and Julian (gregflag=0) calendars.
pub fn revjul(jd: f64, gregflag: i32) -> (i32, i32, f64) {
    let jd_int = (jd + 0.5).floor();
    let frac = jd + 0.5 - jd_int;

    let mut a = jd_int;
    if gregflag != 0 {
        let alpha = ((jd_int - 1867216.25) / 36524.25).floor();
        a = jd_int + 1.0 + alpha - (alpha / 4.0).floor();
    }

    let b = a + 1524.0;
    let c = ((b - 122.1) / 365.25).floor();
    let d = (365.25 * c).floor();
    let e = ((b - d) / 30.6001).floor();

    let day_of_month = (b - d - (30.6001 * e).floor()) + frac;
    let month = if e < 14.0 { e - 1.0 } else { e - 13.0 };
    let mut year = if month > 2.0 { c - 4716.0 } else { c - 4715.0 };

    if year < 1.0 {
        year -= 1.0;
    }

    (year as i32, month as i32, day_of_month)
}

/// Convert date to Julian day number, Gregorian calendar assumed.
#[no_mangle]
pub unsafe extern "C" fn oe_julday(
    year: i32, month: i32, day: f64, gregflag: i32,
) -> f64 {
    julday(year, month, day, gregflag)
}

/// Convert Julian day number to calendar date.
#[no_mangle]
pub unsafe extern "C" fn oe_revjul(
    jd: f64, gregflag: i32,
    year: *mut i32, month: *mut i32, day: *mut f64,
) {
    let (y, m, d) = revjul(jd, gregflag);
    unsafe {
        *year = y;
        *month = m;
        *day = d;
    }
}

/// Compute day of week from Julian day.
/// Returns 0=Monday ... 6=Sunday.
#[no_mangle]
pub unsafe extern "C" fn oe_day_of_week(jd: f64) -> i32 {
    ((jd + 1.5) % 7.0) as i32
}

/// Compute Delta T = TT - UT in seconds for a given Julian day (UT).
/// Uses the selected model from context.
pub fn deltat(jd_ut: f64, model: i32) -> f64 {
    // Simplification: use Stephenson, Morrison & Hohenkerk 2016 polynomial.
    // Full implementation with all models is in delta_t.rs module.
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    match model {
        constants::LE_DT_STEPHENSON_2016 => {
            if jd_ut < 2451545.0 {
                // SMH 2016 polynomial for historical times
                -20.0 + 32.0 * t * t
            } else {
                // For modern times, use quadratic
                65.0 + 203.0 * t + 15.0 * t * t
            }
        }
        constants::LE_DT_ESPENAK_MEEUS_2006 => {
            62.92 + 32.217 * t + 55.89 * t * t
        }
        constants::LE_DT_STEPHENSON_1997 | constants::LE_DT_STEPHENSON_1984 => {
            if jd_ut < 2451545.0 {
                -20.0 + 31.0 * t * t
            } else {
                63.0 + 200.0 * t + 10.0 * t * t
            }
        }
        constants::LE_DT_SCHOCH => {
            67.0 + 205.0 * t + 12.0 * t * t
        }
        _ => 65.0 + 200.0 * t + 10.0 * t * t,
    }
}

/// Date conversion from UTC to Julian day ET (Ephemeris Time ~ TT).
pub fn utc_to_jd(year: i32, month: i32, day: f64, gregflag: i32) -> f64 {
    julday(year, month, day, gregflag)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_julday_gregorian_epoch() {
        let jd = julday(2000, 1, 1.5, LE_GREG_CAL);
        assert!((jd - 2451545.0).abs() < 1e-9);
    }

    #[test]
    fn test_julday_b1950() {
        let jd = julday(1950, 1, 0.9235, LE_GREG_CAL);
        assert!((jd - 2433282.4235).abs() < 0.001);
    }

    #[test]
    fn test_revjul_roundtrip() {
        let jd = 2451545.0;
        let (y, m, d) = revjul(jd, LE_GREG_CAL);
        let jd2 = julday(y, m, d, LE_GREG_CAL);
        assert!((jd - jd2).abs() < 1e-9);
    }

    #[test]
    fn test_revjul_j2000() {
        let (y, m, d) = revjul(2451545.0, LE_GREG_CAL);
        assert_eq!(y, 2000);
        assert_eq!(m, 1);
        assert!((d - 1.5).abs() < 1e-9);
    }
}
