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

/// Fixed star catalog support.
///
/// Provides position calculations for the 200 most astrologically
/// significant stars from the Hipparcos catalog (J2000).
use crate::constants;
use crate::precession::precession_matrix_for_model;

mod data;

/// A fixed star entry in the catalog.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct LeStar {
    /// Bayer/Flamsteed name
    pub name: [i8; 32],
    /// RA J2000 in radians
    pub ra: f64,
    /// Dec J2000 in radians
    pub dec: f64,
    /// Parallax in arcseconds
    pub parallax: f64,
    /// Proper motion in RA (radians/year)
    pub pm_ra: f64,
    /// Proper motion in Dec (radians/year)
    pub pm_dec: f64,
    /// Apparent visual magnitude
    pub mag: f64,
}

/// Number of stars in the embedded catalog.
pub const LE_STAR_COUNT: usize = 200;

/// Get a fixed star position at a given Julian day.
///
/// Applies proper motion and precession from J2000 to the target date.
pub fn star_position(
    star: &LeStar,
    jd_et: f64,
    prec_model: i32,
) -> (f64, f64, f64) {
    let t = (jd_et - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let t_years = t * 100.0;

    // Apply proper motion
    let ra = star.ra + star.pm_ra * t_years;
    let dec = star.dec + star.pm_dec * t_years;

    // Convert to cartesian for precession
    let cos_dec = dec.cos();
    let x = cos_dec * ra.cos();
    let y = cos_dec * ra.sin();
    let z = dec.sin();

    // Apply precession from J2000 to date
    let prec_mat = precession_matrix_for_model(jd_et, prec_model);
    let rotated = prec_mat.transform(&[x, y, z]);

    // Convert back to polar
    let pos = equatorial_cartesian_to_ra_dec_dist(rotated);
    (pos.0, pos.1, pos.2)
}

fn equatorial_cartesian_to_ra_dec_dist(v: [f64; 3]) -> (f64, f64, f64) {
    let dist = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let ra = v[1].atan2(v[0]);
    let dec = (v[2] / dist).asin();
    (ra, dec, dist)
}

/// Normalize star name for comparison: lowercase, strip whitespace.
fn normalize_name(name: &str) -> String {
    name.trim().to_lowercase()
}

/// Look up a star by name (case-insensitive, exact match).
pub fn find_star(name: &str) -> Option<&'static LeStar> {
    let needle = normalize_name(name);
    // Linear search through the catalog
    data::CATALOG.iter().find(|star| {
        let bytes: Vec<u8> = star.name.iter().map(|&b| b as u8).collect();
        let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        let star_name = std::str::from_utf8(&bytes[..end]).unwrap_or("");
        normalize_name(star_name) == needle
    })
}

/// Get star by index (0..LE_STAR_COUNT).
pub fn star_by_index(index: usize) -> Option<&'static LeStar> {
    data::CATALOG.get(index)
}

/// C ABI: get star position at a given Julian day.
#[no_mangle]
pub unsafe extern "C" fn le_fixstar(
    star_name: *const i8,
    jd_et: f64,
    iflag: i32,
    xx: *mut f64,
    serr: *mut i8,
) -> i32 {
    let name = if star_name.is_null() {
        return constants::LE_ERR_INVALID_PARAMS;
    } else {
        let cstr = unsafe { std::ffi::CStr::from_ptr(star_name) };
        match cstr.to_str() {
            Ok(s) => s,
            Err(_) => return constants::LE_ERR_INVALID_PARAMS,
        }
    };

    let star = match find_star(name) {
        Some(s) => s,
        None => {
            if !serr.is_null() { unsafe { *serr = 0; } }
            return constants::ERR_NOT_IMPLEMENTED; // star not found
        }
    };

    let prec_model = if (iflag & constants::LE_FLG_J2000) != 0 {
        constants::LE_PREC_IAU_2006
    } else {
        constants::LE_PREC_IAU_2006
    };

    let (ra, dec, dist) = star_position(star, jd_et, prec_model);

    unsafe {
        *xx.add(0) = ra;
        *xx.add(1) = dec;
        *xx.add(2) = dist;
    }

    if !serr.is_null() { unsafe { *serr = 0; } }
    constants::LE_OK
}

/// C ABI: get star count in catalog.
#[no_mangle]
pub unsafe extern "C" fn le_star_count() -> i32 {
    LE_STAR_COUNT as i32
}

/// C ABI: get star data by index.
#[no_mangle]
pub unsafe extern "C" fn le_star_data(index: i32, star_out: *mut LeStar) -> i32 {
    if index < 0 || index >= LE_STAR_COUNT as i32 {
        return constants::ERR_OUT_OF_RANGE;
    }
    if let Some(star) = star_by_index(index as usize) {
        unsafe { *star_out = *star; }
        constants::LE_OK
    } else {
        constants::ERR_OUT_OF_RANGE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_star_position_no_proper_motion() {
        let star = LeStar {
            name: [0i8; 32],
            ra: 2.530809,
            dec: 1.253107,
            parallax: 0.0075,
            pm_ra: 0.0,
            pm_dec: 0.0,
            mag: 2.0,
        };
        let (ra, dec, _) = star_position(&star, constants::LE_J2000, constants::LE_PREC_IAU_2006);
        // At J2000 with IAU 2006 precession (includes bias), position should be
        // very close to input (within ~2.65 arcsec = ~1.28e-5 rad)
        assert!((ra - 2.530809).abs() < 1e-4);
        assert!((dec - 1.253107).abs() < 1e-4);
    }
}
