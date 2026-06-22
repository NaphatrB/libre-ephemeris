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

/// Solar and lunar eclipse calculations.
///
/// Based on algorithms from Meeus, "Astronomical Algorithms" (2nd ed.),
/// Chapters 54-55. Computes eclipse type and magnitude at a given time.
use crate::constants;
use crate::types::LeVec6;
use crate::analytical;

const AU_KM: f64 = 149_597_870.7;
const SUN_RADIUS_KM: f64 = 695_700.0;
const MOON_RADIUS_KM: f64 = 1737.4;
const EARTH_RADIUS_KM: f64 = 6378.137;

fn get_geocentric_ecliptic(tjd: f64, ipl: i32) -> Result<(f64, f64, f64), i32> {
    let pos = match ipl {
        0 | 1 => {
            if ipl == 0 {
                let earth_helio = analytical::compute_position(tjd, constants::LE_EARTH)?;
                LeVec6::new(
                    -earth_helio.0[0], -earth_helio.0[1], -earth_helio.0[2],
                    -earth_helio.0[3], -earth_helio.0[4], -earth_helio.0[5],
                )
            } else {
                analytical::compute_position(tjd, constants::LE_MOON)?
            }
        }
        _ => return Err(constants::LE_ERR_INVALID_PARAMS),
    };
    let eps0 = (23.439291111111111_f64).to_radians();
    let x = pos.0[0];
    let y = pos.0[1] * eps0.cos() - pos.0[2] * eps0.sin();
    let z = pos.0[1] * eps0.sin() + pos.0[2] * eps0.cos();
    let dist = (x * x + y * y + z * z).sqrt();
    if dist < 1e-15 {
        return Err(constants::ERR_ENGINE);
    }
    let lon = y.atan2(x);
    let lat = (z / dist).asin();
    Ok((lon, lat, dist))
}

fn apparent_radius_rad(distance_au: f64, radius_km: f64) -> f64 {
    let dist_km = distance_au * AU_KM;
    (radius_km / dist_km.max(1.0)).asin()
}

/// Compute solar eclipse parameters at a given time.
///
/// attr array receives:
///   [0] = eclipse type (bitfield: 1=total, 2=annular, 4=partial, 8=central)
///   [1] = magnitude
///   [2] = Sun/Moon apparent diameter ratio
///   [3..5] = Sun geocentric position (lon, lat, dist)
///   [6..8] = Moon geocentric position (lon, lat, dist)
pub fn solar_eclipse_how(tjd: f64, attr: &mut [f64; 20]) -> i32 {
    let sun_geo = match get_geocentric_ecliptic(tjd, 0) {
        Ok(p) => p,
        Err(e) => return e,
    };
    let moon_geo = match get_geocentric_ecliptic(tjd, 1) {
        Ok(p) => p,
        Err(e) => return e,
    };

    let (sun_lon, sun_lat, sun_dist) = sun_geo;
    let (moon_lon, moon_lat, moon_dist) = moon_geo;

    let sep = (moon_lat.sin() * sun_lat.sin() + moon_lat.cos() * sun_lat.cos() * (moon_lon - sun_lon).cos()).acos();

    let sun_rad = apparent_radius_rad(sun_dist, SUN_RADIUS_KM);
    let moon_rad = apparent_radius_rad(moon_dist, MOON_RADIUS_KM);

    let mut etype = 0i32;
    let mag: f64;

    if sep < sun_rad + moon_rad {
        etype |= 4;
        if sep < (sun_rad - moon_rad).abs() {
            if moon_rad > sun_rad {
                etype |= 1;
            } else {
                etype |= 2;
            }
            if moon_rad > sun_rad - sep {
                etype |= 8;
            }
        }
        mag = (sun_rad + moon_rad - sep) / (2.0 * sun_rad);
    } else {
        mag = 0.0;
    }

    attr[0] = etype as f64;
    attr[1] = mag;
    attr[2] = sun_rad / moon_rad;
    attr[3] = sun_lon;
    attr[4] = sun_lat;
    attr[5] = sun_dist;
    attr[6] = moon_lon;
    attr[7] = moon_lat;
    attr[8] = moon_dist;

    0
}

/// Compute lunar eclipse parameters at a given time.
///
/// attr array receives:
///   [0] = eclipse type (bitfield: 1=total, 2=partial, 4=penumbral)
///   [1] = umbral magnitude
///   [2] = penumbral magnitude
///   [3..5] = Moon geocentric position (lon, lat, dist)
///   [6..8] = Earth's shadow center (geocentric ecliptic)
pub fn lunar_eclipse_how(tjd: f64, attr: &mut [f64; 20]) -> i32 {
    let moon_geo = match get_geocentric_ecliptic(tjd, 1) {
        Ok(p) => p,
        Err(e) => return e,
    };
    let sun_geo = match get_geocentric_ecliptic(tjd, 0) {
        Ok(p) => p,
        Err(e) => return e,
    };

    let (moon_lon, moon_lat, moon_dist) = moon_geo;
    let (sun_lon, sun_lat, sun_dist) = sun_geo;

    let shadow_lon = sun_lon + std::f64::consts::PI;
    let shadow_lat = -sun_lat;

    let sep = (moon_lat.sin() * shadow_lat.sin()
        + moon_lat.cos() * shadow_lat.cos() * (moon_lon - shadow_lon).cos())
        .acos();

    let moon_dist_km = moon_dist * AU_KM;
    let sun_dist_km = sun_dist * AU_KM;

    let penumbral_radius = (EARTH_RADIUS_KM + (SUN_RADIUS_KM - EARTH_RADIUS_KM) * (moon_dist_km / sun_dist_km).max(0.0)) / moon_dist_km;
    let umbral_radius = ((EARTH_RADIUS_KM - (SUN_RADIUS_KM - EARTH_RADIUS_KM) * (moon_dist_km / sun_dist_km).max(0.0)).max(0.0)) / moon_dist_km;

    let moon_rad = MOON_RADIUS_KM / moon_dist_km;

    let mut etype = 0i32;
    let umbral_mag: f64;

    if sep < penumbral_radius + moon_rad {
        etype |= 4;
        umbral_mag = (umbral_radius + moon_rad - sep) / (2.0 * moon_rad);
        if umbral_mag > 0.0 {
            etype |= 2;
            if sep + moon_rad < umbral_radius {
                etype |= 1;
            }
        }
    } else {
        umbral_mag = 0.0;
    }

    let penumbral_mag = (penumbral_radius + moon_rad - sep) / (2.0 * moon_rad);

    attr[0] = etype as f64;
    attr[1] = umbral_mag.max(0.0);
    attr[2] = penumbral_mag.max(0.0);
    attr[3] = moon_lon;
    attr[4] = moon_lat;
    attr[5] = moon_dist;
    attr[6] = shadow_lon;
    attr[7] = shadow_lat;
    attr[8] = (EARTH_RADIUS_KM / AU_KM) * (1.0 / moon_dist).max(0.0);

    0
}
