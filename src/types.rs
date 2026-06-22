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

/// 6-vector representing position and velocity in a given coordinate system.
/// The interpretation depends on context: (x,y,z,vx,vy,vz) for cartesian,
/// (lon,lat,dist,lon_dot,lat_dot,dist_dot) for polar, etc.
#[derive(Debug, Clone, Copy, Default)]
#[repr(C)]
pub struct LeVec6(pub [f64; 6]);

impl LeVec6 {
    pub fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self([x, y, z, vx, vy, vz])
    }
    pub fn x(&self) -> f64 { self.0[0] }
    pub fn y(&self) -> f64 { self.0[1] }
    pub fn z(&self) -> f64 { self.0[2] }
    pub fn vx(&self) -> f64 { self.0[3] }
    pub fn vy(&self) -> f64 { self.0[4] }
    pub fn vz(&self) -> f64 { self.0[5] }
}

/// Full ephemeris result: up to 24 doubles covering multiple coordinate systems.
/// Layout: [xyz_eq, xyz_ecl, xyz_eq_of_date, xyz_ecl_of_date]
/// where each is a 6-vector (pos + vel).
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct LePosition {
    pub xx: [f64; 24],
}

impl Default for LePosition {
    fn default() -> Self {
        Self { xx: [0.0; 24] }
    }
}

/// Ephemeris computation flags: bitfield combining OE_FLG_* values
pub type LeFlag = i32;

/// Planet/body identifier (0-17 for fixed bodies, 18+ for asteroids)
pub type LeBody = i32;

/// Julian day number (Ephemeris Time or Terrestrial Time)
pub type LeTime = f64;

/// Geographical latitude in degrees (positive north)
pub type LeLatitude = f64;

/// Geographical longitude in degrees (positive east)
pub type LeLongitude = f64;

/// Altitude above sea level in meters
pub type LeAltitude = f64;

/// Error string buffer (matching C ABI convention)
pub type LeErrBuf = [i8; constants::LE_MAX_ERR_LEN];

/// House system identifier (single character: 'P', 'K', 'E', etc.)
pub type LeHouseSystem = u8;

/// Matrix 3x3 for rotation transformations
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct LeMat3(pub [[f64; 3]; 3]);

impl Default for LeMat3 {
    fn default() -> Self {
        Self([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
}

impl LeMat3 {
    pub fn identity() -> Self { Self::default() }
    pub fn rotate_x(angle_rad: f64) -> Self {
        let (s, c) = angle_rad.sin_cos();
        Self([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])
    }
    pub fn rotate_y(angle_rad: f64) -> Self {
        let (s, c) = angle_rad.sin_cos();
        Self([[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]])
    }
    pub fn rotate_z(angle_rad: f64) -> Self {
        let (s, c) = angle_rad.sin_cos();
        Self([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    }
    pub fn transform(&self, v: &[f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for i in 0..3 {
            out[i] = self.0[i][0] * v[0] + self.0[i][1] * v[1] + self.0[i][2] * v[2];
        }
        out
    }
    pub fn mul(&self, rhs: &LeMat3) -> Self {
        let mut out = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    out[i][j] += self.0[i][k] * rhs.0[k][j];
                }
            }
        }
        Self(out)
    }
    pub fn transpose(&self) -> Self {
        Self([
            [self.0[0][0], self.0[1][0], self.0[2][0]],
            [self.0[0][1], self.0[1][1], self.0[2][1]],
            [self.0[0][2], self.0[1][2], self.0[2][2]],
        ])
    }
}

/// Observer topocentric data
#[derive(Debug, Clone, Copy, Default)]
#[repr(C)]
pub struct LeTopoData {
    pub geolon: f64,
    pub geolat: f64,
    pub geoalt: f64,
}

/// Sidereal (ayanamsa) mode data
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct LeSidData {
    pub mode: i32,
    pub t0: f64,
    pub ayan_t0: f64,
}

impl Default for LeSidData {
    fn default() -> Self {
        Self { mode: constants::LE_SIDM_FAGAN_BRADLEY, t0: 0.0, ayan_t0: 0.0 }
    }
}

/// Ephemeris engine selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LeEngine {
    None,
    Jpl,
    Analytical,
}

/// Physical constants used in calculations
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct LeConst {
    pub speed_of_light: f64,
    pub astronomical_unit: f64,
    pub earth_radius_km: f64,
    pub earth_flattening: f64,
    pub earth_mu: f64,
    pub sun_mu: f64,
    pub j2: f64,
}

impl Default for LeConst {
    fn default() -> Self {
        Self {
            speed_of_light: 173.144633,      // AU/day
            astronomical_unit: 1.495978707e8, // km
            earth_radius_km: 6378.137,
            earth_flattening: 1.0 / 298.257223563,
            earth_mu: 8.997e-10,   // AU^3/day^2 (approx)
            sun_mu: constants::LE_GAUSS_G * constants::LE_GAUSS_G,
            j2: 1.0826359e-3,
        }
    }
}
