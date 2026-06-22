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

/// Analytical ephemeris engine using VSOP87 (planets) and ELP-MPP02 (Moon).
///
/// VSOP87 provides sub-arcsecond positions for all major planets.
/// Time parameter is in millennia from J2000 (TT).
use crate::types::LeVec6;
use crate::constants;

pub mod series;
pub mod planets;
pub mod moon;
pub mod vsop2013;
pub mod elpmpp02;

/// Compute body position using the analytical engine.
///
/// Returns heliocentric (for planets) or geocentric (for Moon) position
/// in equatorial J2000 coordinates, in AU.
pub fn compute_position(jd_et: f64, ipl: i32) -> Result<LeVec6, i32> {
    let t_jcent = (jd_et - constants::LE_J2000) / constants::LE_DAY_PER_CENTURY;

    match ipl {
        2 => Ok(planets::mercury(jd_et)),
        3 => Ok(planets::venus(jd_et)),
        4 => Ok(planets::mars(jd_et)),
        5 => Ok(planets::jupiter(jd_et)),
        6 => Ok(planets::saturn(jd_et)),
        7 => Ok(planets::uranus(jd_et)),
        8 => Ok(planets::neptune(jd_et)),
         9 => Ok(planets::pluto(jd_et)),
         10 => Ok(planets::chiron(jd_et)),
         0 => Ok(LeVec6::new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
        1 => Ok(moon::moon(t_jcent)),
        17 => Ok(planets::earth(jd_et)),
        11 => {
            let earth = planets::earth(jd_et);
            let moon_geo = moon::moon(t_jcent);
            let mu = 0.0123 / (1.0 + 0.0123);
            Ok(LeVec6::new(
                earth.0[0] - moon_geo.0[0] * mu,
                earth.0[1] - moon_geo.0[1] * mu,
                earth.0[2] - moon_geo.0[2] * mu,
                earth.0[3] - moon_geo.0[3] * mu,
                earth.0[4] - moon_geo.0[4] * mu,
                earth.0[5] - moon_geo.0[5] * mu,
            ))
        }
        _ => Err(-1),
    }
}
