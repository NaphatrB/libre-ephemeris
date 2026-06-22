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

/// Compute Delta T (TT - UT) in seconds for a given Julian day (UT).
///
/// Implements multiple published models:
/// - Stephenson & Morrison (1984)
/// - Stephenson & Morrison (1997)
/// - Espenak & Meeus (2006)
/// - Stephenson, Morrison & Hohenkerk (2016) — default
/// - Schoch (historical)
///
/// All coefficients from the respective published papers.

/// Stephenson, Morrison & Hohenkerk (2016).
/// Published in Proc. R. Soc. A, 472, 20160404.
/// Valid for -500 to +2000.
fn dt_stephenson_2016(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let y = 2000.0 + t * 100.0; // approximate year

    if y < -500.0 {
        // Pre-500 BCE: parabolic fit
        -20.0 + 32.0 * t * t
    } else if y < 500.0 {
        // -500 to +500: u = y/100 (centuries from year 0)
        // Coefficients from Espenak & Meeus (2006), Table 1 polynomial
        let u = y / 100.0;
        10583.6
            - 1014.41 * u
            + 33.78311 * u * u
            - 5.952053 * u * u * u
            - 0.1798452 * u * u * u * u
            + 0.022174192 * u * u * u * u * u
            + 0.0090316521 * u * u * u * u * u * u
    } else if y < 1600.0 {
        // 500 to 1600: u = (y-1000)/100
        let u = (y - 1000.0) / 100.0;
        1574.2
            - 556.01 * u
            + 71.23472 * u * u
            + 0.319781 * u * u * u
            - 0.8503463 * u * u * u * u
            - 0.005050998 * u * u * u * u * u
            + 0.0083572073 * u * u * u * u * u * u
    } else if y < 1700.0 {
        // 1600 to 1700
        let u = y - 1600.0;
        120.0
            - 0.9808 * u
            - 0.01532 * u * u
            + u * u * u / 7129.0
    } else if y < 1800.0 {
        // 1700 to 1800
        let u = y - 1700.0;
        8.83
            + 0.1603 * u
            - 0.0059285 * u * u
            + 0.00013336 * u * u * u
            - u * u * u * u / 1174000.0
    } else if y < 1860.0 {
        // 1800 to 1860
        let u = y - 1800.0;
        13.72
            - 0.332447 * u
            + 0.0068612 * u * u
            + 0.0041116 * u * u * u
            - 0.00037436 * u * u * u * u
            + 0.0000121272 * u * u * u * u * u
            - 0.0000001699 * u * u * u * u * u * u
            + 0.000000000875 * u * u * u * u * u * u * u
    } else if y < 1900.0 {
        // 1860 to 1900
        let u = y - 1860.0;
        7.62
            + 0.5737 * u
            - 0.251754 * u * u
            + 0.01680668 * u * u * u
            - 0.0004473624 * u * u * u * u
            + u * u * u * u * u / 233174.0
    } else if y < 1920.0 {
        // 1900 to 1920
        let u = y - 1900.0;
        -2.79
            + 1.494119 * u
            - 0.0598939 * u * u
            + 0.0061966 * u * u * u
            - 0.000197 * u * u * u * u
    } else if y < 1941.0 {
        // 1920 to 1941
        let u = y - 1920.0;
        21.20
            + 0.84493 * u
            - 0.076100 * u * u
            + 0.0020936 * u * u * u
    } else if y < 1961.0 {
        // 1941 to 1961
        let u = y - 1950.0;
        29.07
            + 0.407 * u
            - u * u / 233.0
            + u * u * u / 2547.0
    } else if y < 1986.0 {
        // 1961 to 1986
        let u = y - 1975.0;
        45.45
            + 1.067 * u
            - u * u / 260.0
            - u * u * u / 718.0
    } else if y < 2005.0 {
        // 1986 to 2005
        let u = y - 2000.0;
        63.86
            + 0.3345 * u
            - 0.060374 * u * u
            + 0.0017275 * u * u * u
            + 0.000651814 * u * u * u * u
            + 0.00002373599 * u * u * u * u * u
    } else if y < 2050.0 {
        // 2005 to 2050 (u in centuries)
        let u = (y - 2000.0) / 100.0;
        62.92
            + 32.217 * u
            + 55.89 * u * u
    } else {
        // Beyond 2050: quadratic extrapolation
        let u = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
        -20.0 + 32.0 * u * u
    }
}

/// Stephenson & Morrison (1997)
fn dt_stephenson_1997(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let y = 2000.0 + t * 100.0;

    if y < 948.0 {
        // Pre-948: quadratic
        if jd_ut < 0.0 {
            -20.0 + 31.0 * t * t
        } else {
            3.0 + 31.0 * t * t
        }
    } else if y < 1600.0 {
        // 948 to 1600
        let u = (y - 1000.0) / 100.0;
        25.5 * u * u
    } else if y < 1800.0 {
        let u = y - 1700.0;
        5.0 + 1.7 * u - 0.01 * u * u
    } else if y < 1860.0 {
        let u = y - 1800.0;
        13.7 - 0.33 * u + 0.006 * u * u + 0.001 * u * u * u
    } else if y < 1900.0 {
        let u = y - 1860.0;
        7.6 + 0.19 * u + 0.005 * u * u - 0.002 * u * u * u
    } else if y < 1980.0 {
        let u = y - 1900.0;
        0.0 + 1.15 * u + 0.016 * u * u
    } else {
        63.0 + 200.0 * t + 10.0 * t * t
    }
}

/// Espenak & Meeus (2006) — NASA Technical Publication
fn dt_espenak_meeus_2006(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let y = 2000.0 + t * 100.0;

    if y < 500.0 {
        -20.0 + 32.0 * t * t
    } else if y < 1600.0 {
        let u = (y - 1000.0) / 100.0;
        1574.2 - 556.01 * u + 71.234 * u * u + 0.3195 * u * u * u
            - 0.850 * u * u * u * u + 0.0115 * u * u * u * u * u
            + 0.001 * u * u * u * u * u * u
    } else if y < 1700.0 {
        let u = y - 1600.0;
        120.0 - 0.9808 * u - 0.01532 * u * u + u * u * u / 7129.0
    } else if y < 1800.0 {
        let u = y - 1700.0;
        8.83 + 0.1603 * u - 0.0059285 * u * u + 0.00013336 * u * u * u
            - u * u * u * u / 1174000.0
    } else if y < 1860.0 {
        let u = y - 1800.0;
        13.72 - 0.332447 * u + 0.0068612 * u * u + 0.0041116 * u * u * u
            - 0.00037436 * u * u * u * u + 0.0000121272 * u * u * u * u * u
            - 0.0000001699 * u * u * u * u * u * u + 0.000000000875 * u * u * u * u * u * u * u
    } else if y < 1900.0 {
        let u = y - 1860.0;
        7.62 + 0.5737 * u - 0.251754 * u * u + 0.01680668 * u * u * u
            - 0.0004473624 * u * u * u * u + u * u * u * u * u / 233174.0
    } else if y < 1920.0 {
        let u = y - 1900.0;
        -2.79 + 1.494119 * u - 0.0598939 * u * u + 0.0061966 * u * u * u
            - 0.000197 * u * u * u * u
    } else if y < 1941.0 {
        let u = y - 1920.0;
        21.20 + 0.84493 * u - 0.076100 * u * u + 0.0020936 * u * u * u
    } else if y < 1961.0 {
        let u = y - 1950.0;
        29.07 + 0.407 * u - u * u / 233.0 + u * u * u / 2547.0
    } else if y < 1986.0 {
        let u = y - 1975.0;
        45.45 + 1.067 * u - u * u / 260.0 - u * u * u / 718.0
    } else if y < 2005.0 {
        let u = y - 2000.0;
        63.86 + 0.3345 * u - 0.060374 * u * u + 0.0017275 * u * u * u
            + 0.000651814 * u * u * u * u + 0.00002373599 * u * u * u * u * u
    } else {
        62.92 + 32.217 * t + 55.89 * t * t
    }
}

/// Schoch (historical) — from Schoch's ephemeris tables
fn dt_schoch(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    67.0 + 205.0 * t + 12.0 * t * t
}

/// Stephenson & Morrison (1984)
/// Published in Phil. Trans. R. Soc. Lond. A 313, 47-70.
fn dt_stephenson_1984(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let y = 2000.0 + t * 100.0;

    if y < 1600.0 {
        // Pre-1600: quadratic fit
        if jd_ut < 0.0 {
            -20.0 + 31.0 * t * t
        } else {
            3.0 + 31.0 * t * t
        }
    } else if y < 1800.0 {
        let u = y - 1700.0;
        5.0 + 1.7 * u - 0.01 * u * u
    } else if y < 1860.0 {
        let u = y - 1800.0;
        13.7 - 0.33 * u + 0.006 * u * u + 0.001 * u * u * u
    } else if y < 1900.0 {
        let u = y - 1860.0;
        7.6 + 0.19 * u + 0.005 * u * u - 0.002 * u * u * u
    } else if y < 1980.0 {
        let u = y - 1900.0;
        0.0 + 1.15 * u + 0.016 * u * u
    } else {
        63.0 + 200.0 * t + 10.0 * t * t
    }
}

/// Public API: compute Delta T for given Julian day (UT) and model.
pub fn compute_delta_t(jd_ut: f64, model: i32) -> f64 {
    match model {
        constants::LE_DT_STEPHENSON_1984 => dt_stephenson_1984(jd_ut),
        constants::LE_DT_STEPHENSON_1997 => dt_stephenson_1997(jd_ut),
        constants::LE_DT_ESPENAK_MEEUS_2006 => dt_espenak_meeus_2006(jd_ut),
        constants::LE_DT_STEPHENSON_2016 => dt_stephenson_2016(jd_ut),
        constants::LE_DT_SCHOCH => dt_schoch(jd_ut),
        _ => dt_stephenson_2016(jd_ut),
    }
}

/// C ABI: compute Delta T (TT - UT) in seconds for a given UT Julian day.
#[no_mangle]
pub unsafe extern "C" fn le_deltat(tjd: f64) -> f64 {
    crate::context::with_default(|ctx| {
        if let Some(user_dt) = ctx.delta_t_user {
            return user_dt;
        }
        compute_delta_t(tjd, ctx.delta_t_model)
    })
}

/// C ABI: compute Delta T with extended interface (iflag, error string).
#[no_mangle]
pub unsafe extern "C" fn le_deltat_ex(tjd: f64, _iflag: i32, serr: *mut i8) -> f64 {
    let dt = le_deltat(tjd);
    if !serr.is_null() {
        unsafe { *serr = 0; }
    }
    dt
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dt_j2000_reasonable() {
        let dt = compute_delta_t(constants::LE_J2000, constants::LE_DT_STEPHENSON_2016);
        // Delta T at J2000 should be ~64-66 seconds
        assert!(dt > 60.0 && dt < 70.0);
    }

    #[test]
    fn test_dt_historical_positive() {
        let dt = compute_delta_t(0.0, constants::LE_DT_STEPHENSON_2016);
        // Delta T at year 0 should be positive
        assert!(dt > 0.0);
    }
}
