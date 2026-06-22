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

// Ephemeris computation flag bits — match the original's bit positions
// for C ABI compatibility. These are logical bit assignments, not creative works.
pub const LE_FLG_JPLEPH: i32 = 0x0002;
pub const LE_FLG_SWIEPH: i32 = 0x0004;
pub const LE_FLG_MOSEPH: i32 = 0x0008;
pub const LE_FLG_VSOP2013: i32 = 0x0001;
pub const LE_FLG_HELIO: i32 = 0x0010;
pub const LE_FLG_BARYHEL: i32 = 0x0020;
pub const LE_FLG_TOPOCTR: i32 = 0x0040;
pub const LE_FLG_XYZ: i32 = 0x0080;
pub const LE_FLG_SPEED: i32 = 0x0100;
pub const LE_FLG_NOABERR: i32 = 0x0200;
pub const LE_FLG_NOGDEFL: i32 = 0x0400;
pub const LE_FLG_NOBIRR: i32 = 0x0800;
pub const LE_FLG_J2000: i32 = 0x1000;
pub const LE_FLG_JPLHOR: i32 = 0x2000;
pub const LE_FLG_SIDEREAL: i32 = 0x4000;
pub const LE_FLG_ICRS: i32 = 0x8000;
pub const LE_FLG_EQUATORIAL: i32 = 0x10000;
pub const LE_FLG_ECLIPTIC: i32 = 0x20000;
pub const LE_FLG_TRUE: i32 = 0x40000;
pub const LE_FLG_USER: i32 = 0x80000;
pub const LE_FLG_NONUT: i32 = 0x100000;
pub const LE_FLG_SPEED3: i32 = 0x200000;
pub const LE_FLG_RADIANT: i32 = 0x400000;
pub const LE_FLG_CENTER_BODY: i32 = 0x800000;
pub const LE_FLG_TRUE_VELOCITY: i32 = 0x1000000;

// Planet indices — matching original numbering for C ABI compatibility
pub const LE_SUN: i32 = 0;
pub const LE_MOON: i32 = 1;
pub const LE_MERCURY: i32 = 2;
pub const LE_VENUS: i32 = 3;
pub const LE_MARS: i32 = 4;
pub const LE_JUPITER: i32 = 5;
pub const LE_SATURN: i32 = 6;
pub const LE_URANUS: i32 = 7;
pub const LE_NEPTUNE: i32 = 8;
pub const LE_PLUTO: i32 = 9;
pub const LE_CHIRON: i32 = 10;
pub const LE_MEAN_BARY: i32 = 11;
pub const LE_TRUE_NODE: i32 = 12;
pub const LE_MEAN_NODE: i32 = 13;
pub const LE_MEAN_APOG: i32 = 14;
pub const LE_OSC_APOG: i32 = 15;
pub const LE_INT_APOG: i32 = 16;
pub const LE_EARTH: i32 = 17;

// Fixed planet count
pub const LE_NPLANETS: i32 = 12;
pub const LE_NSIDEREAL_PLANETS: i32 = 19;
pub const LE_NALL_PLANETS: i32 = 19;

// Coordinate system flags for output selection
pub const LE_ECL_GEO: i32 = 256;
pub const LE_ECL_HEL: i32 = 512;
pub const LE_ECL_MAG: i32 = 1024;
pub const LE_ECL_BARY: i32 = 2048;
pub const LE_EQ_GEO: i32 = 4096;
pub const LE_EQ_HEL: i32 = 8192;
pub const LE_HOR: i32 = 16384;

// Sidereal / ayanamsa modes
pub const LE_SIDM_FAGAN_BRADLEY: i32 = 0;
pub const LE_SIDM_LAHIRI: i32 = 1;
pub const LE_SIDM_DELUCE: i32 = 2;
pub const LE_SIDM_RAMAN: i32 = 3;
pub const LE_SIDM_USHASHASHI: i32 = 4;
pub const LE_SIDM_KRISHNAMURTI: i32 = 5;
pub const LE_SIDM_DJWHAL_KHUL: i32 = 6;
pub const LE_SIDM_YUKTESHWAR: i32 = 7;
pub const LE_SIDM_JN_BHASIN: i32 = 8;
pub const LE_SIDM_BABYLONIAN_KUGLER1: i32 = 9;
pub const LE_SIDM_BABYLONIAN_KUGLER2: i32 = 10;
pub const LE_SIDM_BABYLONIAN_KUGLER3: i32 = 11;
pub const LE_SIDM_BABYLONIAN_HUBER: i32 = 12;
pub const LE_SIDM_BABYLONIAN_ETPSC: i32 = 13;
pub const LE_SIDM_ALDEBARAN_15TAU: i32 = 14;
pub const LE_SIDM_HIPPARCHUS: i32 = 15;
pub const LE_SIDM_SASSANIAN: i32 = 16;
pub const LE_SIDM_GALACTIC_CENTER_0SAG: i32 = 17;
pub const LE_SIDM_J2000: i32 = 18;
pub const LE_SIDM_J1900: i32 = 19;
pub const LE_SIDM_B1950: i32 = 20;
pub const LE_SIDM_SURYASIDDHANTA: i32 = 21;
pub const LE_SIDM_ARYABHATA: i32 = 22;
pub const LE_SIDM_SS_CHITRAPAKSHA: i32 = 23;
pub const LE_SIDM_SS_REVATI: i32 = 24;
pub const LE_SIDM_TRUE_CITRA: i32 = 25;
pub const LE_SIDM_TRUE_REVATI: i32 = 26;
pub const LE_SIDM_TRUE_PUSHYA: i32 = 27;
pub const LE_SIDM_TRUE_ASHVINI: i32 = 28;
pub const LE_SIDM_TRUE_MAGHA: i32 = 29;
pub const LE_SIDM_TRUE_MULA: i32 = 30;
pub const LE_SIDM_GALCENT_CENTER_0SAG: i32 = 31;
pub const LE_SIDM_GALCENT_MULA: i32 = 32;
pub const LE_SIDM_GALCENT_RADIO_0SAG: i32 = 33;
pub const LE_SIDM_GALCENT_0CAP: i32 = 34;
pub const LE_SIDM_TRUE_ANTARES: i32 = 35;
pub const LE_SIDM_TRUE_FOMALHAUT: i32 = 36;
pub const LE_SIDM_VALENTINE: i32 = 37;
pub const LE_SIDM_USER: i32 = 38;
pub const LE_SIDM_LARGE: i32 = 39;
pub const LE_SIDM_TRUE_ALDEBARAN: i32 = 40;
pub const LE_SIDM_STEINER_SET: i32 = 41;
pub const LE_SIDM_RAMAN_A: i32 = 42;
pub const LE_SIDM_SURYAPAKSHA_A: i32 = 43;
pub const LE_SIDM_TRUE_VEGA: i32 = 44;
pub const LE_SIDM_TRUE_ZETA_PSC: i32 = 45;
pub const LE_SIDM_SURYAPAKSHA_B: i32 = 46;
pub const LE_NMODES_AYANAMSA: i32 = 47;

// Precession model IDs
pub const LE_PREC_IAU_1976: i32 = 0;
pub const LE_PREC_LASKAR_1986: i32 = 1;
pub const LE_PREC_WILLIAMS_1994: i32 = 2;
pub const LE_PREC_SIMON_1994: i32 = 3;
pub const LE_PREC_IAU_2000: i32 = 4;
pub const LE_PREC_BRETAGNON_2003: i32 = 5;
pub const LE_PREC_IAU_2006: i32 = 6;
pub const LE_PREC_VONDRAK_2011: i32 = 7;
pub const LE_PREC_OWEN_1990: i32 = 8;
pub const LE_PREC_NEWCOMB: i32 = 9;
pub const LE_PREC_IAU_2006_REDUCED: i32 = 10;
pub const LE_NMODES_PREC: i32 = 11;

// Nutation model IDs
pub const LE_NUT_IAU_1980: i32 = 0;
pub const LE_NUT_IAU_1980_HERRING: i32 = 1;
pub const LE_NUT_IAU_2000A: i32 = 2;
pub const LE_NUT_IAU_2000B: i32 = 3;
pub const LE_NUT_WOOLARD_1953: i32 = 4;
pub const LE_NMODES_NUT: i32 = 5;

// Delta T model IDs
pub const LE_DT_STEPHENSON_1984: i32 = 0;
pub const LE_DT_STEPHENSON_1997: i32 = 1;
pub const LE_DT_ESPENAK_MEEUS_2006: i32 = 2;
pub const LE_DT_STEPHENSON_2016: i32 = 3;
pub const LE_DT_SCHOCH: i32 = 4;
pub const LE_DT_USER: i32 = 5;
pub const LE_NMODES_DT: i32 = 6;

// Sidereal time model IDs
pub const LE_ST_IAU_1976: i32 = 0;
pub const LE_ST_IAU_2006: i32 = 1;
pub const LE_ST_IERS_2010: i32 = 2;
pub const LE_ST_LONG_TERM: i32 = 3;
pub const LE_NMODES_ST: i32 = 4;

// Frame bias model IDs
pub const LE_BIAS_NONE: i32 = 0;
pub const LE_BIAS_IAU_2000: i32 = 1;
pub const LE_BIAS_IAU_2006: i32 = 2;
pub const LE_NMODES_BIAS: i32 = 3;

// JPL Horizons approximation mode
pub const LE_HOR_APPROX_NONE: i32 = 0;
pub const LE_HOR_APPROX_STANDARD: i32 = 1;
pub const LE_HOR_APPROX_REFINED: i32 = 2;

// Error codes (negative = error, 0+ = success with info)
pub const LE_OK: i32 = 0;
pub const LE_ERR: i32 = -1;
pub const LE_ERR_INVALID_PARAMS: i32 = -2;
pub const ERR_FILE_NOT_FOUND: i32 = -3;
pub const ERR_OUT_OF_RANGE: i32 = -4;
pub const ERR_NO_EPHEMERIS: i32 = -5;
pub const ERR_IO: i32 = -6;
pub const ERR_ENGINE: i32 = -7;
pub const ERR_INVALID_PLANET: i32 = -8;
pub const ERR_INVALID_FLAG: i32 = -9;
pub const ERR_MEMORY: i32 = -10;
pub const ERR_NOT_IMPLEMENTED: i32 = -11;
pub const ERR_UNKNOWN_BODY: i32 = -12;

// Output array indices for le_calc etc.
// The output xx array is 6 doubles for position+velocity.
// With LE_SIDEREAL or additional flags, up to 24 doubles.
pub const LE_X: usize = 0;
pub const LE_Y: usize = 1;
pub const LE_Z: usize = 2;
pub const LE_VX: usize = 3;
pub const LE_VY: usize = 4;
pub const LE_VZ: usize = 5;

pub const LE_RA: usize = 0;
pub const LE_DEC: usize = 1;
pub const LE_DIST: usize = 2;
pub const LE_RA_DOT: usize = 3;
pub const LE_DEC_DOT: usize = 4;
pub const LE_DIST_DOT: usize = 5;

pub const LE_LON: usize = 0;
pub const LE_LAT: usize = 1;
pub const LE_DISTANCE: usize = 2;
pub const LE_LON_DOT: usize = 3;
pub const LE_LAT_DOT: usize = 4;
pub const LE_DISTANCE_DOT: usize = 5;

// Max string length for error messages
pub const LE_MAX_ERR_LEN: usize = 256;

// Standard gravitational constant (in AU^3 / day^2)
pub const LE_GAUSS_G: f64 = 0.01720209895;

// Julian day constants
pub const LE_J1970: f64 = 2440587.5;
pub const LE_J2000: f64 = 2451545.0;
pub const LE_B1950: f64 = 2433282.4235;
pub const LE_J1900: f64 = 2415020.0;

// Days per unit time
pub const LE_DAY_PER_YEAR: f64 = 365.25;
pub const LE_DAY_PER_CENTURY: f64 = 36525.0;
pub const LE_DAY_PER_MILLENNIUM: f64 = 365250.0;

// Angles
pub const LE_DEG: f64 = 0.017453292519943295; // π/180
pub const LE_RAD: f64 = 57.29577951308232;    // 180/π
pub const LE_ARCSEC: f64 = 4.84813681109536e-6; // π/(180*3600)
pub const LE_ARCMIN: f64 = 2.908882086657216e-4; // π/(180*60)

// Speed of light in AU/day
pub const LE_CLIGHT: f64 = 173.1446326846693;
