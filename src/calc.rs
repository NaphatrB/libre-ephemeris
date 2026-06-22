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

/// Core ephemeris calculation engine.
///
/// Orchestrates ephemeris engine selection, coordinate transformations,
/// corrections (precession, nutation, aberration, light deflection, bias),
/// and frame conversions.
use crate::constants;
use crate::types::*;
use crate::context;

fn resolve_body_index(ipl: i32) -> i32 {
    ipl // pass-through for standard bodies
}

/// Main ephemeris calculation.
///
/// For a given body at a given Julian date, compute the 6-vector (position + velocity)
/// in the requested coordinate system, applying all requested corrections.
///
/// Returns up to 24 doubles in `xx`:
///   xx[0..5]   = position/velocity in selected output system
///   Sidereal mode adds 6 more doubles in xx[6..11]
///   etc.
pub fn calc(
    tjd: f64,
    ipl: i32,
    iflag: i32,
    xx: &mut [f64; 24],
    serr: &mut String,
) -> i32
{
    context::with_default(|ctx| {
        xx.fill(0.0);

        let body = resolve_body_index(ipl);
        if body < 0 || body > 17 {
            serr.push_str("invalid planet index");
            return constants::ERR_INVALID_PLANET;
        }

        // Determine engine
        let use_jpl = (iflag & constants::LE_FLG_JPLEPH) != 0;
        let use_vsop2013 = (iflag & constants::LE_FLG_VSOP2013) != 0;
        let _use_analytical = !use_jpl;

        // Compute raw heliocentric J2000 equatorial position
        let raw_pos = if use_jpl {
            #[cfg(feature = "jpl")]
            {
                crate::jpl::compute_position(tjd, body)
            }
            #[cfg(not(feature = "jpl"))]
            {
                Err(constants::ERR_ENGINE)
            }
        } else if use_vsop2013 {
            #[cfg(feature = "analytical")]
            {
                use crate::analytical::vsop2013;
                // Moon via ELP-MPP02 if available
                if body == constants::LE_MOON {
                    if let Some(ref reader) = ctx.elpmpp02_reader {
                        reader.compute_moon(tjd).map_err(|_| constants::ERR_ENGINE)
                    } else {
                        crate::analytical::compute_position(tjd, body)
                    }
                } else {
                    let body_map = |b: i32| -> Option<usize> {
                        match b {
                            2 => Some(1),  // Mercury
                            3 => Some(2),  // Venus
                            4 => Some(4),  // Mars
                            5 => Some(5),  // Jupiter
                            6 => Some(6),  // Saturn
                            7 => Some(7),  // Uranus
                            8 => Some(8),  // Neptune
                            9 => Some(9),  // Pluto
                            11 | 17 => Some(3), // Earth-Moon barycenter / Earth
                            _ => None,
                        }
                    };
                    if let Some(vsop_idx) = body_map(body) {
                        if let Some(ref reader) = ctx.vsop2013_readers[vsop_idx - 1] {
                            vsop2013::compute_position_vsop2013(tjd, body, reader)
                        } else {
                            crate::analytical::compute_position(tjd, body)
                        }
                    } else {
                        Err(constants::ERR_ENGINE)
                    }
                }
            }
            #[cfg(not(feature = "analytical"))]
            {
                Err(constants::ERR_ENGINE)
            }
        } else {
            #[cfg(feature = "analytical")]
            {
                crate::analytical::compute_position(tjd, body)
            }
            #[cfg(not(feature = "analytical"))]
            {
                Err(constants::ERR_ENGINE)
            }
        };

        let raw = match raw_pos {
            Ok(p) => p,
            Err(e) => {
                serr.push_str("ephemeris computation failed");
                return e;
            }
        };

        // Get Earth position (for barycentric correction and aberration).
        // Needed for any body that undergoes frame conversion or aberration,
        // or for Moon when converting to/from barycentric/heliocentric.
        let need_earth = body == constants::LE_MOON
            || (body != constants::LE_MEAN_BARY && body != constants::LE_EARTH);
        let earth_pos = if need_earth {
            Some(get_earth_position(tjd, use_jpl, ctx))
        } else {
            None
        };

        // Start with raw heliocentric equatorial J2000 position
        let pos = [raw.0[0], raw.0[1], raw.0[2]];
        let vel = [raw.0[3], raw.0[4], raw.0[5]];

        // Determine output frame
        let output_helio = (iflag & constants::LE_FLG_HELIO) != 0;
        let output_bary = (iflag & constants::LE_FLG_BARYHEL) != 0;
        // Default (neither): geocentric

        // For planets, corrections work in heliocentric frame.
        // For Moon, keep geocentric (analytical engine already returns this).
        // We apply frame conversion at the end.

        // Store intermediate J2000 equatorial position
        for i in 0..6 {
            xx[i] = if i < 3 { pos[i] } else { vel[i - 3] };
        }

        // Apply corrections unless explicitly disabled
        let mut corrected_pos = pos;
        let mut corrected_vel = vel;

        // Frame bias (ICRS -> FK5 or vice versa)
        if (iflag & constants::LE_FLG_NOBIRR) == 0 {
            corrected_pos = crate::bias::apply_frame_bias(&corrected_pos, ctx.bias_model);
            corrected_vel = crate::bias::apply_frame_bias(&corrected_vel, ctx.bias_model);
        }

        // Precession (J2000 -> date).
        // When LE_FLG_J2000 is set, skip precession to keep output in J2000 frame.
        if (iflag & constants::LE_FLG_J2000) == 0 {
            let prec_mat = crate::precession::precession_matrix_for_model(tjd, ctx.prec_model);
            corrected_pos = prec_mat.transform(&corrected_pos);
            corrected_vel = prec_mat.transform(&corrected_vel);
        }

        // Nutation (true equator of date)
        if (iflag & constants::LE_FLG_NONUT) == 0 {
            let nut = crate::nutation::compute_nutation(tjd, ctx.nut_model);
            let nut_mat = crate::nutation::nutation_matrix(&nut);
            corrected_pos = nut_mat.transform(&corrected_pos);
            corrected_vel = nut_mat.transform(&corrected_vel);
        }

        // Gravitational light deflection by Sun
        if (iflag & constants::LE_FLG_NOGDEFL) == 0 && body != constants::LE_SUN {
            if let Some(ep) = &earth_pos {
                let defl = crate::deflection::light_deflection(
                    &corrected_pos,
                    &[ep.0[0], ep.0[1], ep.0[2]],
                    &[0.0; 3], // Sun at barycenter
                );
                corrected_pos[0] += defl[0];
                corrected_pos[1] += defl[1];
                corrected_pos[2] += defl[2];
            }
        }

        // Annual aberration
        if (iflag & constants::LE_FLG_NOABERR) == 0 {
            if let Some(ep) = &earth_pos {
                let dist = (corrected_pos[0] * corrected_pos[0]
                    + corrected_pos[1] * corrected_pos[1]
                    + corrected_pos[2] * corrected_pos[2])
                    .sqrt();
                if dist > 1e-10 {
                    let dir = [
                        corrected_pos[0] / dist,
                        corrected_pos[1] / dist,
                        corrected_pos[2] / dist,
                    ];
                    let earth_vel = [ep.0[3], ep.0[4], ep.0[5]];
                    let ab_corr = crate::aberration::aberration_correction(&dir, &earth_vel);
                    corrected_pos[0] += ab_corr[0] * dist;
                    corrected_pos[1] += ab_corr[1] * dist;
                    corrected_pos[2] += ab_corr[2] * dist;
                }
            }
        }

        // Topocentric correction
        if (iflag & constants::LE_FLG_TOPOCTR) != 0 && ctx.topo_set {
            // Compute local sidereal time
            let lst = crate::houses::ramc_from_lst(
                greenwich_sidereal_time(tjd) + ctx.topo.geolon
            );
            let obs_geoc = crate::topocentric::geodetic_to_geocentric(
                ctx.topo.geolon,
                ctx.topo.geolat,
                ctx.topo.geoalt,
                lst * constants::LE_DEG,
            );
            // Get Earth position
            if let Some(ep) = &earth_pos {
                let body_earth = [
                    corrected_pos[0] - ep.0[0],
                    corrected_pos[1] - ep.0[1],
                    corrected_pos[2] - ep.0[2],
                ];
                let topo = crate::topocentric::topocentric_correction(&body_earth, &obs_geoc);
                corrected_pos = [
                    topo[0] + ep.0[0],
                    topo[1] + ep.0[1],
                    topo[2] + ep.0[2],
                ];
            }
        }

        // Sidereal (ayanamsa) correction
        if (iflag & constants::LE_FLG_SIDEREAL) != 0 {
            let aya = crate::ayanamsa::compute_ayanamsa(tjd, &ctx.sid_mode) * constants::LE_DEG;
            // Rotate around ecliptic pole (z-axis of ecliptic)
            // Convert to ecliptic, subtract ayanamsa, convert back
            let eps0 = crate::transform::mean_obliquity_iau2006(0.0);
            let ecl = crate::transform::equatorial_to_ecliptic_cartesian(
                corrected_pos[0], corrected_pos[1], corrected_pos[2], eps0
            );
            let ecl_vel = crate::transform::equatorial_to_ecliptic_cartesian(
                corrected_vel[0], corrected_vel[1], corrected_vel[2], eps0
            );
            let (lon, lat, dist) = crate::transform::ecliptic_cartesian_to_polar(ecl[0], ecl[1], ecl[2]);
            let (_lon_dot, _lat_dot, _dist_dot) = crate::transform::ecliptic_cartesian_to_polar(ecl_vel[0], ecl_vel[1], ecl_vel[2]);
            let lon_sid = lon - aya;
            let (sin_lon, cos_lon) = lon_sid.sin_cos();
            let cos_lat = lat.cos();
            let sid_pos = [dist * cos_lat * cos_lon, dist * cos_lat * sin_lon, dist * lat.sin()];
            corrected_pos = crate::transform::ecliptic_to_equatorial_cartesian(
                sid_pos[0], sid_pos[1], sid_pos[2], eps0
            );
        }

        // Convert to requested output frame FIRST, then apply coordinate system conversion.
        //
        // Raw positions from engine:
        //   - planets: heliocentric
        //   - Moon: geocentric
        //   - Sun: heliocentric (origin)
        //
        // Conversions:
        //   HELIO: keep heliocentric (for Moon: geo → helio by subtracting Earth)
        //   BARYHEL: add Earth position → approximate barycentric
        //   Default (geocentric): subtract Earth → position from Earth (for Moon: keep geo)
        let frame_pos: [f64; 3] = if output_bary {
            if body == constants::LE_MOON {
                if let Some(ep) = &earth_pos {
                    [corrected_pos[0] + ep.0[0], corrected_pos[1] + ep.0[1], corrected_pos[2] + ep.0[2]]
                } else { corrected_pos }
            } else if let Some(ep) = &earth_pos {
                [corrected_pos[0] + ep.0[0], corrected_pos[1] + ep.0[1], corrected_pos[2] + ep.0[2]]
            } else { corrected_pos }
        } else if output_helio {
            if body == constants::LE_MOON {
                if let Some(ep) = &earth_pos {
                    [corrected_pos[0] - ep.0[0], corrected_pos[1] - ep.0[1], corrected_pos[2] - ep.0[2]]
                } else { corrected_pos }
            } else { corrected_pos }
        } else {
            // Default: geocentric
            if body != constants::LE_MOON {
                if let Some(ep) = &earth_pos {
                    [corrected_pos[0] - ep.0[0], corrected_pos[1] - ep.0[1], corrected_pos[2] - ep.0[2]]
                } else { corrected_pos }
            } else { corrected_pos }
        };

        // Write position and velocity to output
        for i in 0..6 {
            xx[i] = if i < 3 { frame_pos[i] } else { corrected_vel[i - 3] };
        }

        // Coordinate system conversions (applied on top of frame-converted position)
        if (iflag & constants::LE_FLG_ECLIPTIC) != 0 {
            // Use J2000 obliquity when output is in J2000 frame, else use date's obliquity
            let eps = if (iflag & constants::LE_FLG_J2000) != 0 {
                crate::transform::mean_obliquity_iau2006(0.0)
            } else {
                crate::transform::mean_obliquity_iau2006(
                    (tjd - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY
                )
            };
            let ecl_pos = crate::transform::equatorial_to_ecliptic_cartesian(
                xx[0], xx[1], xx[2], eps
            );
            let (lon, lat, dist) = crate::transform::ecliptic_cartesian_to_polar(ecl_pos[0], ecl_pos[1], ecl_pos[2]);
            if (iflag & constants::LE_FLG_XYZ) != 0 {
                xx[0] = ecl_pos[0]; xx[1] = ecl_pos[1]; xx[2] = ecl_pos[2];
            } else {
                xx[0] = lon * constants::LE_RAD;
                xx[1] = lat * constants::LE_RAD;
                xx[2] = dist;
            }
        } else if (iflag & constants::LE_FLG_EQUATORIAL) != 0 || (iflag & constants::LE_FLG_XYZ) == 0 {
            // Polar equatorial output (in degrees)
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            let ra = xx[1].atan2(xx[0]) * constants::LE_RAD;
            let dec = (xx[2] / dist).asin() * constants::LE_RAD;
            xx[0] = ra; xx[1] = dec; xx[2] = dist;
        }

        constants::LE_OK
    })
}

fn get_earth_position(tjd: f64, use_jpl: bool, ctx: &context::LeContext) -> LeVec6 {
    // Compute Earth position for use as barycentric offset
    if use_jpl {
        #[cfg(feature = "jpl")]
        {
            crate::jpl::compute_position(tjd, constants::LE_EARTH).unwrap_or_default()
        }
        #[cfg(not(feature = "jpl"))]
        { LeVec6::default() }
    } else if ctx.vsop2013_readers.len() > 2 && ctx.vsop2013_readers[2].is_some() {
        #[cfg(feature = "analytical")]
        {
            let reader = ctx.vsop2013_readers[2].as_ref().unwrap();
            crate::analytical::vsop2013::compute_position_vsop2013(tjd, constants::LE_EARTH, reader)
                .unwrap_or_else(|_| crate::analytical::compute_position(tjd, constants::LE_EARTH).unwrap_or_default())
        }
        #[cfg(not(feature = "analytical"))]
        { LeVec6::default() }
    } else {
        #[cfg(feature = "analytical")]
        {
            crate::analytical::compute_position(tjd, constants::LE_EARTH).unwrap_or_default()
        }
        #[cfg(not(feature = "analytical"))]
        { LeVec6::default() }
    }
}

fn greenwich_sidereal_time(jd_ut: f64) -> f64 {
    let t = (jd_ut - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;
    let t2 = t * t;
    let t3 = t2 * t;
    let gmst_sec = 24110.54841 + 8640184.812866 * t + 0.093104 * t2 - 0.0000062 * t3;
    let gmst_deg = (gmst_sec / 240.0) % 360.0;
    if gmst_deg < 0.0 { gmst_deg + 360.0 } else { gmst_deg }
}

/// C ABI: main ephemeris calculation (ET input).
#[no_mangle]
pub unsafe extern "C" fn le_calc(
    tjd: f64,
    ipl: i32,
    iflag: i32,
    xx: *mut f64,
    serr: *mut i8,
) -> i32 {
    let mut xx_arr = [0.0_f64; 24];
    let mut error = String::new();

    let rc = calc(tjd, ipl, iflag, &mut xx_arr, &mut error);

    for i in 0..24 {
        unsafe { *xx.add(i) = xx_arr[i]; }
    }

    // Copy error string to C buffer
    if !serr.is_null() && !error.is_empty() {
        let bytes = error.as_bytes();
        let len = bytes.len().min(constants::LE_MAX_ERR_LEN - 1);
        unsafe {
            std::ptr::copy_nonoverlapping(bytes.as_ptr() as *const i8, serr, len);
            *serr.add(len) = 0;
        }
    } else if !serr.is_null() {
        unsafe { *serr = 0; }
    }

    rc
}

/// C ABI: main ephemeris calculation (UT input).
///
/// Computes the position and velocity of a celestial body at a given
/// Universal Time. Internally converts UT to ET using Delta T.
///
/// # Arguments
/// * `tjd_ut` - Julian day in Universal Time
/// * `ipl` - Planet index (0=Sun, 1=Moon, 2=Mercury, ..., 9=Pluto, 10=Chiron, 17=Earth)
/// * `iflag` - Bitmask of OE_FLG_* flags controlling output frame, coordinate system, corrections
/// * `xx` - Output array of 24 doubles: [x, y, z, vx, vy, vz, ...]
/// * `serr` - Error string buffer (256 bytes), zeroed on success
///
/// # Returns
/// 0 on success, negative on error.
///
/// # Flags
/// - `LE_FLG_XYZ`: Cartesian output (default: polar)
/// - `LE_FLG_J2000`: J2000 frame (default: date)
/// - `LE_FLG_ECLIPTIC`: Ecliptic coordinates (default: equatorial)
/// - `LE_FLG_HELIO`: Heliocentric frame (default: geocentric)
/// - `LE_FLG_BARYHEL`: Barycentric frame
/// - `LE_FLG_NOABERR`: Skip annual aberration correction
/// - `LE_FLG_NOGDEFL`: Skip gravitational deflection
/// - `LE_FLG_NOBIRR`: Skip frame bias
/// - `LE_FLG_NONUT`: Skip nutation
/// - `LE_FLG_SPEED`: Include velocity in output
/// - `LE_FLG_SIDEREAL`: Apply sidereal (ayanamsa) correction
/// - `LE_FLG_TOPOCTR`: Topocentric position
/// - `LE_FLG_VSOP2013`: Use VSOP2013/ELP-MPP02 file-based engine
/// - `LE_FLG_JPLEPH`: Use JPL DE ephemeris
#[no_mangle]
pub unsafe extern "C" fn le_calc_ut(
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
    xx: *mut f64,
    serr: *mut i8,
) -> i32 {
    // Convert UT to ET by adding Delta T
    let dt_seconds = crate::delta_t::le_deltat(tjd_ut);
    let tjd_et = tjd_ut + dt_seconds / 86400.0;
    le_calc(tjd_et, ipl, iflag, xx, serr)
}

/// C ABI: solar eclipse calculation.
#[no_mangle]
pub unsafe extern "C" fn oe_sol_eclipse_how(
    tjd: f64,
    _iflag: i32,
    attr: *mut f64,
    serr: *mut i8,
) -> i32 {
    let mut attr_arr = [0.0_f64; 20];
    let rc = crate::eclipse::solar_eclipse_how(tjd, &mut attr_arr);
    for i in 0..20 {
        unsafe { *attr.add(i) = attr_arr[i]; }
    }
    if !serr.is_null() { unsafe { *serr = 0; } }
    rc
}

/// C ABI: lunar eclipse calculation.
#[no_mangle]
pub unsafe extern "C" fn oe_lun_eclipse_how(
    tjd: f64,
    _iflag: i32,
    attr: *mut f64,
    serr: *mut i8,
) -> i32 {
    let mut attr_arr = [0.0_f64; 20];
    let rc = crate::eclipse::lunar_eclipse_how(tjd, &mut attr_arr);
    for i in 0..20 {
        unsafe { *attr.add(i) = attr_arr[i]; }
    }
    if !serr.is_null() { unsafe { *serr = 0; } }
    rc
}
