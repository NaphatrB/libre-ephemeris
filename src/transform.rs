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

/// Convert equatorial polar coordinates (RA, Dec, distance in radians/AU)
/// to equatorial cartesian (x, y, z in AU).
pub fn equatorial_polar_to_cartesian(ra: f64, dec: f64, dist: f64) -> [f64; 3] {
    let d = dist * dec.cos();
    [d * ra.cos(), d * ra.sin(), dist * dec.sin()]
}

/// Convert equatorial cartesian to equatorial polar.
pub fn equatorial_cartesian_to_polar(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let dist = (x * x + y * y + z * z).sqrt();
    let ra = y.atan2(x);
    let dec = (z / dist).asin();
    (ra, dec, dist)
}

/// Mean obliquity of the ecliptic (IAU 1976), in radians.
/// Lieske et al. (1977), A&A 58, 1-16.
/// ε = 23°26'21.448" − 46.8150" T − 0.00059" T² + 0.001813" T³
/// where T is in Julian centuries from J2000.
pub fn mean_obliquity_iau1976(t_jcent: f64) -> f64 {
    let eps0 = 23.43929111111111; // degrees at J2000
    let dpsi = -46.8150 * t_jcent - 0.00059 * t_jcent * t_jcent + 0.001813 * t_jcent * t_jcent * t_jcent;
    (eps0 + dpsi / 3600.0) * constants::LE_DEG
}

/// Mean obliquity of the ecliptic (IAU 2006), in radians.
/// Capitaine et al. (2003), A&A 412, 567-586.
/// ε = 84381.406″ − 46.836769″ T − 0.0001831″ T² + 0.00200340″ T³
///     − 0.576×10⁻⁶ T⁴ − 4.34×10⁻⁸ T⁵
/// where T is in Julian centuries from J2000.
pub fn mean_obliquity_iau2006(t_jcent: f64) -> f64 {
    let eps0 = 84381.406; // arcseconds at J2000
    let dpsi = -46.836769 * t_jcent
        - 0.0001831 * t_jcent * t_jcent
        + 0.00200340 * t_jcent * t_jcent * t_jcent
        - 0.576e-6 * t_jcent * t_jcent * t_jcent * t_jcent
        - 4.34e-8 * t_jcent * t_jcent * t_jcent * t_jcent * t_jcent;
    (eps0 + dpsi) * (std::f64::consts::PI / (180.0 * 3600.0))
}

/// Convert equatorial cartesian to ecliptic cartesian.
/// Rotation by obliquity of the ecliptic around the x-axis.
pub fn equatorial_to_ecliptic_cartesian(
    x: f64, y: f64, z: f64, eps: f64,
) -> [f64; 3] {
    let cos_eps = eps.cos();
    let sin_eps = eps.sin();
    let xe = x;
    let ye = y * cos_eps + z * sin_eps;
    let ze = -y * sin_eps + z * cos_eps;
    [xe, ye, ze]
}

/// Convert ecliptic cartesian to equatorial cartesian.
pub fn ecliptic_to_equatorial_cartesian(
    x: f64, y: f64, z: f64, eps: f64,
) -> [f64; 3] {
    let cos_eps = eps.cos();
    let sin_eps = eps.sin();
    let xq = x;
    let yq = y * cos_eps - z * sin_eps;
    let zq = y * sin_eps + z * cos_eps;
    [xq, yq, zq]
}

/// Convert ecliptic polar (lon, lat in radians, distance in AU)
/// to ecliptic cartesian.
pub fn ecliptic_polar_to_cartesian(lon: f64, lat: f64, dist: f64) -> [f64; 3] {
    let d = dist * lat.cos();
    [d * lon.cos(), d * lon.sin(), dist * lat.sin()]
}

/// Convert ecliptic cartesian to ecliptic polar.
pub fn ecliptic_cartesian_to_polar(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let dist = (x * x + y * y + z * z).sqrt();
    let lon = y.atan2(x);
    let lat = (z / dist).asin();
    (lon, lat, dist)
}

/// Equatorial to horizontal coordinate transformation.
/// phi = geographical latitude, theta = local sidereal time (both in radians).
pub fn equatorial_to_horizontal(
    ra: f64, dec: f64, phi: f64, theta: f64,
) -> (f64, f64) {
    let h = theta - ra; // hour angle
    let alt = (phi.sin() * dec.sin() + phi.cos() * dec.cos() * h.cos()).asin();
    let az = -(dec.sin() - phi.sin() * alt).atan2(dec.cos() * h.sin());
    // Normalize azimuth: 0 = north, 90 = east
    let az = (az + std::f64::consts::PI) % (2.0 * std::f64::consts::PI);
    (az, alt)
}

/// Horizontal to equatorial coordinate transformation.
pub fn horizontal_to_equatorial(
    az: f64, alt: f64, phi: f64, theta: f64,
) -> (f64, f64) {
    let h = (alt.sin() * phi.sin() - az.cos() * alt.cos() * phi.cos()).atan2(
        -(az.sin() * alt.cos()),
    );
    let dec = (alt.sin() * phi.cos() + az.cos() * alt.cos() * phi.sin()).asin();
    let ra = theta - h;
    (ra, dec)
}

/// Split decimal degrees into degrees, minutes, seconds.
pub fn split_deg(ddeg: f64) -> (i32, i32, f64) {
    let deg = ddeg.trunc() as i32;
    let rem = (ddeg - deg as f64).abs() * 60.0;
    let min = rem.trunc() as i32;
    let sec = (rem - min as f64) * 60.0;
    (deg, min, sec)
}

/// Normalize degrees to [0, 360).
pub fn csnorm(d: f64) -> f64 {
    let d = d % 360.0;
    if d < 0.0 { d + 360.0 } else { d }
}

/// Difference of two angles in degrees (result in [0, 360)).
pub fn difdeg2n(d1: f64, d2: f64) -> f64 {
    let d = (d1 - d2) % 360.0;
    if d < -180.0 { d + 360.0 } else if d >= 180.0 { d - 360.0 } else { d }
}

/// C ABI: coordinate transformation between systems.
/// iflag bits select source and target.
/// 0: equatorial J2000 -> ecliptic J2000
/// 1: ecliptic J2000 -> equatorial J2000
/// 2: equatorial J2000 -> equatorial of date
/// etc.
#[no_mangle]
pub unsafe extern "C" fn le_cotrans(
    x: *const f64, y: *const f64, z: *const f64,
    eps: f64,
    xout: *mut f64, yout: *mut f64, zout: *mut f64,
) {
    let (xi, yi, zi) = unsafe { (*x, *y, *z) };
    let result = if (desired_flags() & 1) != 0 {
        let r = equatorial_to_ecliptic_cartesian(xi, yi, zi, eps);
        (r[0], r[1], r[2])
    } else {
        let r = ecliptic_to_equatorial_cartesian(xi, yi, zi, eps);
        (r[0], r[1], r[2])
    };
    unsafe {
        *xout = result.0;
        *yout = result.1;
        *zout = result.2;
    }
}

// Placeholder: retrieve desired conversion direction.
// In the full implementation this would be driven by le_cotrans_sp or
// a dedicated flag parameter.
fn desired_flags() -> i32 { 0 }

/// C ABI: coordinate transformation for speed vectors (6-vector).
#[no_mangle]
pub unsafe extern "C" fn le_cotrans_sp(
    x: *const f64, y: *const f64, z: *const f64,
    eps: f64,
    xout: *mut f64, yout: *mut f64, zout: *mut f64,
) {
    le_cotrans(x, y, z, eps, xout, yout, zout);
}

/// C ABI: split decimal degrees.
#[no_mangle]
pub unsafe extern "C" fn le_split_deg(
    ddeg: f64, _roundflag: i32,
    ideg: *mut i32, imin: *mut i32, isec: *mut f64, dsecfrac: *mut f64,
) {
    let (d, m, s) = split_deg(ddeg);
    unsafe {
        *ideg = d;
        *imin = m;
        *isec = s;
        if !dsecfrac.is_null() {
            *dsecfrac = s.fract();
        }
    }
}

/// C ABI: normalize angle to [0, 360).
#[no_mangle]
pub unsafe extern "C" fn le_csnorm(d: f64) -> f64 {
    csnorm(d)
}

/// C ABI: difference of two angles.
#[no_mangle]
pub unsafe extern "C" fn le_difdeg2n(d1: f64, d2: f64) -> f64 {
    difdeg2n(d1, d2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_equatorial_to_ecliptic_roundtrip() {
        let eps = mean_obliquity_iau2006(0.0);
        let x = 1.0;
        let y = 0.0;
        let z = 0.0;
        let ecl = equatorial_to_ecliptic_cartesian(x, y, z, eps);
        let eq = ecliptic_to_equatorial_cartesian(ecl[0], ecl[1], ecl[2], eps);
        assert!((eq[0] - x).abs() < 1e-15);
        assert!((eq[1] - y).abs() < 1e-15);
        assert!((eq[2] - z).abs() < 1e-15);
    }

    #[test]
    fn test_split_deg() {
        let (d, m, s) = split_deg(23.5);
        assert_eq!(d, 23);
        assert_eq!(m, 30);
        assert!((s - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_norm() {
        assert!((csnorm(370.0) - 10.0).abs() < 1e-12);
        assert!((csnorm(-10.0) - 350.0).abs() < 1e-12);
    }
}
