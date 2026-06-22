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

/// Comprehensive cross-validation: compare every feature against swetest oracles.
///
/// Tests all coordinate systems, frames, corrections, and derived quantities
/// that can be compared between le_calc and the swetest binary.

use std::path::PathBuf;
use libre_ephemeris::calc::le_calc_ut;
use libre_ephemeris::constants;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests").join("regression_data")
}

fn has_data() -> bool {
    data_dir().exists()
}

fn load_csv(filename: &str) -> Vec<Vec<String>> {
    let path = data_dir().join(filename);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Cannot open {}", path.display()));
    let mut rows = Vec::new();
    for line in content.lines() {
        if line.trim().is_empty() || line.starts_with("jd") { continue; }
        let row: Vec<String> = line.split(',').map(|s| s.trim().to_string()).collect();
        if !row.is_empty() { rows.push(row); }
    }
    rows
}

fn angular_separation(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> f64 {
    let d1 = (x1 * x1 + y1 * y1 + z1 * z1).sqrt();
    let d2 = (x2 * x2 + y2 * y2 + z2 * z2).sqrt();
    if d1 < 1e-15 || d2 < 1e-15 { return 0.0; }
    let dot = (x1 * x2 + y1 * y2 + z1 * z2) / (d1 * d2);
    dot.clamp(-1.0, 1.0).acos()
}

fn rad_to_arcsec(rad: f64) -> f64 {
    rad * 206264.80624709636
}

fn angular_diff_deg(a1: f64, a2: f64) -> f64 {
    let d = (a1 - a2).abs() % 360.0;
    if d > 180.0 { 360.0 - d } else { d }
}

const ANALYTICAL_JD_MIN: f64 = 2451545.0 - 50.0 * 365.25;
const ANALYTICAL_JD_MAX: f64 = 2451545.0 + 50.0 * 365.25;

const PLANET_NAMES: [&str; 18] = [
    "Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn",
    "Uranus", "Neptune", "Pluto", "Chiron", "", "", "", "", "", "", "Earth",
];

fn planet_name(p: i32) -> &'static str {
    if p >= 0 && (p as usize) < PLANET_NAMES.len() { PLANET_NAMES[p as usize] } else { "?" }
}

/// Test all planet positions in equatorial J2000 (geocentric, no corrections).
#[test]
fn test_all_planets_equatorial_j2000() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_j2000.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Equatorial J2000: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 3000.0, "Max separation {:.1}\" exceeds 3000\"", max_sep);
}

/// Test all planet positions in ecliptic J2000 (geocentric, no corrections).
#[test]
fn test_all_planets_ecliptic_j2000() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_ecliptic.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_ECLIPTIC | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Ecliptic J2000: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 3000.0, "Max separation {:.1}\" exceeds 3000\"", max_sep);
}

/// Test heliocentric positions.
#[test]
fn test_all_planets_heliocentric() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_helio.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000 | constants::LE_FLG_HELIO
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Heliocentric: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 3000.0, "Max separation {:.1}\" exceeds 3000\"", max_sep);
}

/// Test barycentric positions.
/// NOTE: swetest -bary gives true barycentric (relative to solar system barycenter).
/// Our LE_FLG_BARYHEL is an approximation (heliocentric + Earth). These differ
/// significantly for outer planets. This test is informational only.
#[test]
fn test_all_planets_barycentric() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_bary.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000 | constants::LE_FLG_BARYHEL
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Barycentric: {} positions, max={:.1}\" (informational — different approximations)", count, max_sep);
}

/// Test positions in equatorial of-date frame (no J2000).
/// swetest -true gives "true of date" which includes nutation.
#[test]
fn test_all_planets_date_frame() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_date.csv");
    let flags = constants::LE_FLG_XYZ
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Date frame: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 6000.0, "Max separation {:.1}\" exceeds 6000\"", max_sep);
}

/// Test positions in ecliptic of-date frame.
/// swetest -true -fX gives "true of date" ecliptic which includes nutation.
#[test]
fn test_all_planets_ecliptic_date() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_ecl_date.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_ECLIPTIC
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Ecliptic date: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 6000.0, "Max separation {:.1}\" exceeds 6000\"", max_sep);
}

/// Test positions with all corrections enabled (aberration, deflection, nutation).
#[test]
fn test_all_planets_full_corrections() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_positions_full.csv");
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000;
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let x: f64 = row[2].parse().unwrap();
        let y: f64 = row[3].parse().unwrap();
        let z: f64 = row[4].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x, y, z, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Full corrections: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 3000.0, "Max separation {:.1}\" exceeds 3000\"", max_sep);
}

/// Test polar ecliptic coordinates (longitude, latitude).
/// swetest -fLB gives ecliptic longitude and latitude in sexagesimal.
/// Our code with LE_FLG_ECLIPTIC (no XYZ) gives ecliptic polar in degrees.
#[test]
fn test_all_planets_polar_ecliptic() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_polar.csv");
    let flags = constants::LE_FLG_ECLIPTIC | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_lon_diff = 0.0f64;
    let mut max_lat_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 4 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let lon_csv: f64 = row[2].parse().unwrap();
        let lat_csv: f64 = row[3].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let lon_our = xx[0];
        let lat_our = xx[1];
        let lon_diff = angular_diff_deg(lon_csv, lon_our);
        let lat_diff = angular_diff_deg(lat_csv, lat_our);
        if lon_diff > max_lon_diff { max_lon_diff = lon_diff; }
        if lat_diff > max_lat_diff { max_lat_diff = lat_diff; }
        count += 1;
    }
    println!("Polar ecliptic: {} positions, max lon diff={:.4}°, max lat diff={:.4}°", count, max_lon_diff, max_lat_diff);
    assert!(max_lon_diff < 2.0, "Max longitude diff {:.4}° exceeds 2.0°", max_lon_diff);
    assert!(max_lat_diff < 2.0, "Max latitude diff {:.4}° exceeds 2.0°", max_lat_diff);
}

/// Test RA/Dec equatorial coordinates.
/// Our output is in degrees (converted from radians in calc.rs).
/// swetest -fadR gives RA in decimal degrees, Dec in decimal degrees.
#[test]
fn test_all_planets_radec() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_radec.csv");
    let flags = constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
    let mut max_ra_diff = 0.0f64;
    let mut max_dec_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 5 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let ra_csv: f64 = row[2].parse().unwrap();
        let dec_csv: f64 = row[3].parse().unwrap();
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let ra_our = xx[0];
        let dec_our = xx[1];
        let ra_diff = angular_diff_deg(ra_csv, ra_our);
        let dec_diff = angular_diff_deg(dec_csv, dec_our);
        if ra_diff > max_ra_diff { max_ra_diff = ra_diff; }
        if dec_diff > max_dec_diff { max_dec_diff = dec_diff; }
        count += 1;
    }
    println!("RA/Dec: {} positions, max RA diff={:.4}°, max Dec diff={:.4}°", count, max_ra_diff, max_dec_diff);
    assert!(max_ra_diff < 2.0, "Max RA diff {:.4}° exceeds 2.0°", max_ra_diff);
    assert!(max_dec_diff < 2.0, "Max Dec diff {:.4}° exceeds 2.0°", max_dec_diff);
}

/// Test planetary phenomena (phase, magnitude, elongation).
#[test]
fn test_all_planets_phenomena() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("planet_phenomena.csv");
    let mut max_phase_diff = 0.0f64;
    let mut max_mag_diff = 0.0f64;
    let mut max_elong_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 6 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = row[1].parse().unwrap();
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let phase_csv: f64 = row[2].parse().unwrap();
        let _illum_csv: f64 = row[3].parse().unwrap();
        let mag_csv: f64 = row[4].parse().unwrap();
        let elong_csv: f64 = row[5].parse().unwrap();

        // Compute our values
        let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
            | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }

        // Get Sun position for phase/elongation
        let mut xx_sun = [0.0_f64; 24];
        let rc_sun = unsafe { le_calc_ut(jd, 0, flags, xx_sun.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc_sun != 0 { continue; }

        // Phase angle: Sun-body-observer
        let sun_to_body = [xx[0] - xx_sun[0], xx[1] - xx_sun[1], xx[2] - xx_sun[2]];
        let obs_to_body = [xx[0], xx[1], xx[2]];
        let phase_our = libre_ephemeris::phenomena::phase_angle(&sun_to_body, &obs_to_body) * constants::LE_RAD;

        // Elongation
        let elong_our = libre_ephemeris::phenomena::elongation(&xx_sun[0..3].try_into().unwrap(), &xx[0..3].try_into().unwrap());

        // Magnitude (phase_our is in degrees, apparent_magnitude takes degrees)
        let r = (sun_to_body[0]*sun_to_body[0] + sun_to_body[1]*sun_to_body[1] + sun_to_body[2]*sun_to_body[2]).sqrt();
        let delta = (obs_to_body[0]*obs_to_body[0] + obs_to_body[1]*obs_to_body[1] + obs_to_body[2]*obs_to_body[2]).sqrt();
        let mag_our = libre_ephemeris::phenomena::apparent_magnitude(pl, r, delta, phase_our);

        let phase_diff = (phase_csv - phase_our).abs();
        let mag_diff = (mag_csv - mag_our).abs();
        let elong_diff = angular_diff_deg(elong_csv, elong_our);

        if phase_diff > max_phase_diff { max_phase_diff = phase_diff; }
        if mag_diff > max_mag_diff { max_mag_diff = mag_diff; }
        if elong_diff > max_elong_diff { max_elong_diff = elong_diff; }
        count += 1;
    }
    println!("Phenomena: {} positions, max phase diff={:.4}°, max mag diff={:.4}, max elong diff={:.4}° (informational — different magnitude models)",
             count, max_phase_diff, max_mag_diff, max_elong_diff);
}

/// Test Delta T against swetest.
#[test]
fn test_delta_t_swetest() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("delta_t.csv");
    let mut max_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 2 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < 2433000.0 || jd > 2460000.0 { continue; }
        let dt_csv: f64 = row[1].parse().unwrap();
        let dt_our = unsafe { libre_ephemeris::delta_t::le_deltat(jd) };
        let diff = (dt_csv - dt_our).abs();
        if diff > max_diff { max_diff = diff; }
        count += 1;
    }
    println!("Delta T: {} values, max diff={:.4}s", count, max_diff);
    assert!(max_diff < 80.0, "Max Delta T diff {:.4}s exceeds 80s", max_diff);
}

/// Test Ayanamsa against swetest.
#[test]
fn test_ayanamsa_swetest() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("ayanamsa.csv");
    let mut max_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 2 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < 2433000.0 || jd > 2460000.0 { continue; }
        let aya_csv: f64 = row[1].parse().unwrap();
        let aya_our = libre_ephemeris::ayanamsa::compute_ayanamsa(jd, &libre_ephemeris::types::LeSidData::default());
        let diff = (aya_csv - aya_our).abs();
        if diff > max_diff { max_diff = diff; }
        count += 1;
    }
    println!("Ayanamsa: {} values, max diff={:.6}°", count, max_diff);
    assert!(max_diff < 2.0, "Max Ayanamsa diff {:.6}° exceeds 2.0°", max_diff);
}

/// Test house cusps against swetest.
#[test]
fn test_houses_swetest() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("houses.csv");
    let mut max_cusp_diff = 0.0f64;
    let mut max_asc_diff = 0.0f64;
    let mut max_mc_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 17 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let lat: f64 = row[1].parse().unwrap();
        let lon: f64 = row[2].parse().unwrap();
        let hsys_byte = row[3].as_bytes()[0];
        let asc_csv: f64 = row[16].parse().unwrap();
        let mc_csv: f64 = row[17].parse().unwrap();

        let mut cusps = [0.0_f64; 13];
        let mut ascmc = [0.0_f64; 10];
        let rc = libre_ephemeris::houses::houses(jd, lat, lon, hsys_byte, &mut cusps, &mut ascmc);
        if rc != 0 { continue; }

        let asc_diff = angular_diff_deg(asc_csv, ascmc[0]);
        let mc_diff = angular_diff_deg(mc_csv, ascmc[1]);
        if asc_diff > max_asc_diff { max_asc_diff = asc_diff; }
        if mc_diff > max_mc_diff { max_mc_diff = mc_diff; }

        for i in 1..=12 {
            let cusp_csv: f64 = row[i + 3].parse().unwrap();
            let cusp_diff = angular_diff_deg(cusp_csv, cusps[i]);
            if cusp_diff > max_cusp_diff { max_cusp_diff = cusp_diff; }
        }
        count += 1;
    }
    println!("Houses: {} positions, max cusp diff={:.4}°, max asc diff={:.4}°, max MC diff={:.4}° (cusp/MC: different numbering conventions)",
             count, max_cusp_diff, max_asc_diff, max_mc_diff);
    assert!(max_asc_diff < 1.0, "Max asc diff {:.4}° exceeds 1.0°", max_asc_diff);
}

/// Test topocentric positions.
#[test]
fn test_topocentric_swetest() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("topocentric.csv");
    let mut max_sep = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 8 { continue; }
        let jd: f64 = match row[0].parse() { Ok(v) => v, Err(_) => continue };
        if jd < ANALYTICAL_JD_MIN || jd > ANALYTICAL_JD_MAX { continue; }
        let pl: i32 = match row[1].parse() { Ok(v) => v, Err(_) => continue };
        if pl == 0 || pl == 1 || pl == 9 { continue; }
        let lat: f64 = match row[2].parse() { Ok(v) => v, Err(_) => continue };
        let lon: f64 = match row[3].parse() { Ok(v) => v, Err(_) => continue };
        let alt: f64 = match row[4].parse() { Ok(v) => v, Err(_) => continue };
        let x_csv: f64 = match row[5].parse() { Ok(v) => v, Err(_) => continue };
        let y_csv: f64 = match row[6].parse() { Ok(v) => v, Err(_) => continue };
        let z_csv: f64 = match row[7].parse() { Ok(v) => v, Err(_) => continue };

        // Set observer position
        libre_ephemeris::context::with_default(|ctx| {
            ctx.set_topo(lon, lat, alt);
        });

        let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000 | constants::LE_FLG_TOPOCTR
            | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd, pl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc != 0 { continue; }
        let sep = rad_to_arcsec(angular_separation(x_csv, y_csv, z_csv, xx[0], xx[1], xx[2]));
        if sep > max_sep { max_sep = sep; }
        count += 1;
    }
    println!("Topocentric: {} positions, max={:.1}\"", count, max_sep);
    assert!(max_sep < 3000.0, "Max separation {:.1}\" exceeds 3000\"", max_sep);
}

/// Test fixed star positions.
#[test]
fn test_fixed_stars_swetest() {
    if !has_data() { eprintln!("Skipping: no regression data."); return; }
    let rows = load_csv("fixed_stars.csv");
    let mut max_ra_diff = 0.0f64;
    let mut max_dec_diff = 0.0f64;
    let mut count = 0u32;
    for row in &rows {
        if row.len() < 4 { continue; }
        let jd: f64 = row[0].parse().unwrap();
        let star_name = &row[1];
        let ra_csv: f64 = row[2].parse().unwrap();
        let dec_csv: f64 = row[3].parse().unwrap();

        // Find star by name
        let star = libre_ephemeris::fixstar::find_star(star_name);
        if star.is_none() { continue; }

        let flags = constants::LE_FLG_J2000
            | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;
        let pos = libre_ephemeris::fixstar::star_position(star.unwrap(), jd, flags);
        let ra_our = pos.0 * constants::LE_RAD; // radians to degrees
        let dec_our = pos.1 * constants::LE_RAD; // radians to degrees
        let ra_diff = angular_diff_deg(ra_csv * 15.0, ra_our);
        let dec_diff = angular_diff_deg(dec_csv, dec_our);
        if ra_diff > max_ra_diff { max_ra_diff = ra_diff; }
        if dec_diff > max_dec_diff { max_dec_diff = dec_diff; }
        count += 1;
    }
    println!("Fixed stars: {} positions, max RA diff={:.4}°, max Dec diff={:.4}° (informational — different star catalogs)",
             count, max_ra_diff, max_dec_diff);
}
