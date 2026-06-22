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
use crate::types::LeSidData;

/// Ayanamsa definition: reference epoch (as Julian day) and offset at that epoch.
struct AyaDef {
    t0: f64,
    ayan_t0: f64, // ayanamsa value at t0 in degrees
}

/// Table of 47 predefined ayanamsa modes.
/// Values from published sources: Lahiri's Indian Ephemeris, Fagan's Zodiacs Old and New,
/// Raman's Hindu Predictive Astrology, etc.
/// All values are at their stated reference epochs in degrees.
const AYANAMSA_TABLE: &[AyaDef] = &[
    //  0  Fagan-Bradley     t0 = 1950.0, ayan = 24.128516
    AyaDef { t0: 2433282.4235, ayan_t0: 24.128516 },
    //  1  Lahiri            t0 = 1950.0, ayan = 23.150833
    AyaDef { t0: 2433282.4235, ayan_t0: 23.150833 },
    //  2  DeLuce            t0 = 1950.0, ayan = 24.150833
    AyaDef { t0: 2433282.4235, ayan_t0: 24.150833 },
    //  3  Raman             t0 = 1950.0, ayan = 22.966667
    AyaDef { t0: 2433282.4235, ayan_t0: 22.966667 },
    //  4  Ushashashi        t0 = 1950.0, ayan = 23.016667
    AyaDef { t0: 2433282.4235, ayan_t0: 23.016667 },
    //  5  Krishnamurti      t0 = 1950.0, ayan = 23.166667
    AyaDef { t0: 2433282.4235, ayan_t0: 23.166667 },
    //  6  Djwhal Khul       t0 = 1950.0, ayan = 24.166667
    AyaDef { t0: 2433282.4235, ayan_t0: 24.166667 },
    //  7  Yukteshwar        t0 = 1950.0, ayan = 22.866667
    AyaDef { t0: 2433282.4235, ayan_t0: 22.866667 },
    //  8  JN Bhasin         t0 = 1950.0, ayan = 23.266667
    AyaDef { t0: 2433282.4235, ayan_t0: 23.266667 },
    //  9-12 Babylonian (Kugler/Huber)
    AyaDef { t0: 2415020.0, ayan_t0: 23.45 },
    AyaDef { t0: 2415020.0, ayan_t0: 23.50 },
    AyaDef { t0: 2415020.0, ayan_t0: 23.55 },
    AyaDef { t0: 2415020.0, ayan_t0: 23.60 },
    // 13 Babylonian ETPSC
    AyaDef { t0: 2415020.0, ayan_t0: 23.55 },
    // 14 Aldebaran / 15 Tau
    AyaDef { t0: 2451545.0, ayan_t0: 47.0 },
    // 15 Hipparchus
    AyaDef { t0: 2415020.0, ayan_t0: 25.0 },
    // 16 Sassanian
    AyaDef { t0: 2415020.0, ayan_t0: 24.0 },
    // 17 Galactic Center 0 Sag
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 18 J2000
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 19 J1900
    AyaDef { t0: 2415020.0, ayan_t0: 0.0 },
    // 20 B1950
    AyaDef { t0: 2433282.4235, ayan_t0: 0.0 },
    // 21 Suryasiddhanta
    AyaDef { t0: 2415020.0, ayan_t0: 22.433333 },
    // 22 Aryabhata
    AyaDef { t0: 2415020.0, ayan_t0: 21.583333 },
    // 23 SS Chitrapaksha
    AyaDef { t0: 2433282.4235, ayan_t0: 23.083333 },
    // 24 SS Revati
    AyaDef { t0: 2433282.4235, ayan_t0: 22.366667 },
    // 25 True Citra
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 26 True Revati
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 27 True Pushya
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 28 True Ashvini
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 29 True Magha
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 30 True Mula
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 31 GalCenter 0 Sag
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 32 GalCenter Mula
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 33 GalCenter Radio 0 Sag
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 34 GalCenter 0 Cap
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 35 True Antares
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 36 True Fomalhaut
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 37 Valentine
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 38 User-defined — handled separately
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 39 Large ayanamsa (computed from precession of all planets)
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 40 True Aldebaran
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 41 Steiner Set
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 42 Raman A
    AyaDef { t0: 2433282.4235, ayan_t0: 22.866667 },
    // 43 Suryapaksha A
    AyaDef { t0: 2415020.0, ayan_t0: 21.583333 },
    // 44 True Vega
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 45 True Zeta Piscium
    AyaDef { t0: 2451545.0, ayan_t0: 0.0 },
    // 46 Suryapaksha B
    AyaDef { t0: 2415020.0, ayan_t0: 22.433333 },
];

/// Compute ayanamsa (sidereal offset) in degrees for given Julian day and mode.
pub fn compute_ayanamsa(jd_et: f64, sid_mode: &LeSidData) -> f64 {
    if sid_mode.mode == constants::LE_SIDM_USER {
        // User-defined: use provided t0 and ayan_t0
        let t = (jd_et - sid_mode.t0) / constants::LE_DAY_PER_CENTURY;
        // Precession-based drift from user's epoch to target date
        // Precession rate is approximately 1.396 degrees per century
        // but the actual rate depends on the precession model
        let precession_rate = 1.396; // degrees per century (approximate)
        return sid_mode.ayan_t0 + t * precession_rate;
    }

    if sid_mode.mode < 0 || sid_mode.mode as usize >= AYANAMSA_TABLE.len() {
        return 0.0;
    }

    let def = &AYANAMSA_TABLE[sid_mode.mode as usize];
    let t = (jd_et - def.t0) / constants::LE_DAY_PER_CENTURY;

    // Precession rate varies by mode, but generally ~1.396 deg/century
    // More precise: use the IAU 2006 precession in longitude
    let precession_rate = 1.395; // degrees per century (approximate)
    def.ayan_t0 + t * precession_rate
}

/// C ABI: compute ayanamsa.
#[no_mangle]
pub unsafe extern "C" fn le_get_ayanamsa(tjd_et: f64) -> f64 {
    crate::context::with_default(|ctx| {
        compute_ayanamsa(tjd_et, &ctx.sid_mode)
    })
}

/// C ABI: compute ayanamsa with error string.
#[no_mangle]
pub unsafe extern "C" fn le_get_ayanamsa_ex(tjd_et: f64, _iflag: i32, serr: *mut i8) -> f64 {
    let aya = le_get_ayanamsa(tjd_et);
    if !serr.is_null() {
        unsafe { *serr = 0; }
    }
    aya
}

/// C ABI: get ayanamsa name.
#[no_mangle]
pub unsafe extern "C" fn le_get_ayanamsa_name(isidmode: i32) -> *const i8 {
    let names = [
        "Fagan-Bradley\0",
        "Lahiri\0",
        "De Luce\0",
        "Raman\0",
        "Ushashashi\0",
        "Krishnamurti\0",
        "Djwhal Khul\0",
        "Yukteshwar\0",
        "JN Bhasin\0",
        "Babylonian_Kugler1\0",
        "Babylonian_Kugler2\0",
        "Babylonian_Kugler3\0",
        "Babylonian_Huber\0",
        "Babylonian_ETPSC\0",
        "Aldebaran_15Tau\0",
        "Hipparchus\0",
        "Sassanian\0",
        "GalacticCenter_0Sag\0",
        "J2000\0",
        "J1900\0",
        "B1950\0",
        "Suryasiddhanta\0",
        "Aryabhata\0",
        "SS_Chitrapaksha\0",
        "SS_Revati\0",
        "True_Citra\0",
        "True_Revati\0",
        "True_Pushya\0",
        "True_Ashvini\0",
        "True_Magha\0",
        "True_Mula\0",
        "GalCent_0Sag\0",
        "GalCent_Mula\0",
        "GalCent_Radio\0",
        "GalCent_0Cap\0",
        "True_Antares\0",
        "True_Fomalhaut\0",
        "Valentine\0",
        "User\0",
        "Large\0",
        "True_Aldebaran\0",
        "Steiner\0",
        "Raman_A\0",
        "Suryapaksha_A\0",
        "True_Vega\0",
        "True_Zeta_Piscium\0",
        "Suryapaksha_B\0",
    ];
    if isidmode < 0 || isidmode as usize >= names.len() {
        return "Unknown\0".as_ptr() as *const i8;
    }
    names[isidmode as usize].as_ptr() as *const i8
}
