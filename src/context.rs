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

/// Thread-local context manager for the ephemeris library.
///
/// Manages global state that would otherwise be module-level: ephemeris file paths,
/// observer location (topocentric), model selections (precession, nutation, sidereal,
/// frame bias), and user-defined Delta T overrides. Provides both a Rust API
/// (`with_default`, `read_default`, `LeContext` methods) and C ABI wrappers
/// (`le_set_ephe_path`, `le_set_topo`, etc.).
#[cfg(feature = "analytical")]
use crate::analytical::vsop2013::Vsop2013Reader;
#[cfg(feature = "analytical")]
use crate::analytical::elpmpp02::ElpMpp02Reader;
#[cfg(feature = "analytical")]
use crate::analytical::elpmpp02::ElpParas;
use crate::constants;
use crate::types::*;
use std::cell::RefCell;

/// Internal ephemeris computation context.
/// Manages state that would otherwise be global: ephemeris paths,
/// open file handles, observer position, model selection, etc.
#[derive(Debug)]
pub struct LeContext {
    pub ephe_paths: Vec<String>,
    pub jpl_filename: Option<String>,
    pub jpl_file_is_open: bool,
    #[cfg(feature = "analytical")]
    pub vsop2013_readers: Vec<Option<Vsop2013Reader>>,
    #[cfg(feature = "analytical")]
    pub elpmpp02_reader: Option<ElpMpp02Reader>,
    #[cfg(feature = "analytical")]
    pub elpmpp02_paras: ElpParas,
    pub topo: LeTopoData,
    pub topo_set: bool,
    pub sid_mode: LeSidData,
    pub delta_t_user: Option<f64>,
    pub delta_t_model: i32,
    pub tide_acc: f64,
    pub prec_model: i32,
    pub nut_model: i32,
    pub st_model: i32,
    pub bias_model: i32,
    pub hor_approx: i32,
    pub physical_constants: LeConst,
}

impl LeContext {
    /// Create a new context with default settings.
    pub fn new() -> Self {
        Self {
            ephe_paths: Vec::new(),
            jpl_filename: None,
            jpl_file_is_open: false,
            #[cfg(feature = "analytical")]
            vsop2013_readers: (0..9).map(|_| None).collect(),
            #[cfg(feature = "analytical")]
            elpmpp02_reader: None,
            #[cfg(feature = "analytical")]
            elpmpp02_paras: ElpParas::default(),
            topo: LeTopoData::default(),
            topo_set: false,
            sid_mode: LeSidData::default(),
            delta_t_user: None,
            delta_t_model: constants::LE_DT_STEPHENSON_2016,
            tide_acc: 0.0,
            prec_model: constants::LE_PREC_IAU_2006,
            nut_model: constants::LE_NUT_IAU_2000B,
            st_model: constants::LE_ST_IAU_2006,
            bias_model: constants::LE_BIAS_IAU_2006,
            hor_approx: constants::LE_HOR_APPROX_NONE,
            physical_constants: LeConst::default(),
        }
    }

    /// Set ephemeris search paths (colon/semicolon-separated).
    pub fn set_ephe_path(&mut self, path: &str) {
        self.ephe_paths.clear();
        for p in path.split(&[':', ';'][..]) {
            let trimmed = p.trim();
            if !trimmed.is_empty() {
                self.ephe_paths.push(trimmed.to_string());
            }
        }
    }

    /// Set the JPL DE ephemeris binary file path.
    pub fn set_jpl_file(&mut self, fname: &str) {
        self.jpl_filename = Some(fname.to_string());
    }

    /// Set observer geocentric position (longitude deg, latitude deg, altitude m).
    pub fn set_topo(&mut self, lon: f64, lat: f64, alt: f64) {
        self.topo = LeTopoData { geolon: lon, geolat: lat, geoalt: alt };
        self.topo_set = true;
    }

    /// Set sidereal (ayanamsa) mode.
    pub fn set_sid_mode(&mut self, mode: i32, t0: f64, ayan_t0: f64) {
        self.sid_mode = LeSidData { mode, t0, ayan_t0 };
    }

    /// Override Delta T (TT - UT) in seconds. None = use model.
    pub fn set_delta_t_user(&mut self, dt: f64) {
        self.delta_t_user = Some(dt);
    }

    /// Set tidal acceleration parameter (arcsec/cy²).
    pub fn set_tid_acc(&mut self, acc: f64) {
        self.tide_acc = acc;
    }

    /// Load a VSOP2013 data file for a given planet (1-9).
    #[cfg(feature = "analytical")]
    pub fn set_vsop2013_file(&mut self, planet: usize, path: &str) -> Result<(), String> {
        if planet < 1 || planet > 9 {
            return Err("Planet must be 1-9".to_string());
        }
        let reader = Vsop2013Reader::open(path)?;
        self.vsop2013_readers[planet - 1] = Some(reader);
        Ok(())
    }

    /// Load ELP-MPP02 data files from a directory.
    #[cfg(feature = "analytical")]
    pub fn set_elpmpp02_dir(&mut self, dir: &str) -> Result<(), String> {
        let reader = ElpMpp02Reader::open_with_paras(dir, &self.elpmpp02_paras)?;
        self.elpmpp02_reader = Some(reader);
        Ok(())
    }

    /// Set ELP-MPP02 adjustable parameters (for LLR vs DE405 fit, etc.).
    /// Takes effect on the next call to set_elpmpp02_dir().
    #[cfg(feature = "analytical")]
    pub fn set_elpmpp02_paras(&mut self, paras: &ElpParas) {
        self.elpmpp02_paras = *paras;
    }

    /// Select precession, nutation, sidereal time, and frame bias models.
    pub fn set_astro_models(&mut self, prec: i32, nut: i32, st: i32, bias: i32) {
        self.prec_model = prec;
        self.nut_model = nut;
        self.st_model = st;
        self.bias_model = bias;
    }
}

impl Default for LeContext {
    fn default() -> Self {
        Self::new()
    }
}

// Thread-local default context
thread_local! {
    static DEFAULT_CONTEXT: RefCell<LeContext> = RefCell::new(LeContext::new());
}

/// Access the thread-local default context.
pub fn with_default<F, T>(f: F) -> T
where
    F: FnOnce(&mut LeContext) -> T,
{
    DEFAULT_CONTEXT.with(|ctx| f(&mut ctx.borrow_mut()))
}

/// Read the thread-local default context.
pub fn read_default<F, T>(f: F) -> T
where
    F: FnOnce(&LeContext) -> T,
{
    DEFAULT_CONTEXT.with(|ctx| f(&ctx.borrow()))
}

/// C ABI: set ephemeris search paths (colon/semicolon-separated).
#[no_mangle]
pub unsafe extern "C" fn le_set_ephe_path(path: *const i8) {
    if path.is_null() { return; }
    let s = unsafe { std::ffi::CStr::from_ptr(path) }.to_str().unwrap_or("");
    with_default(|ctx| ctx.set_ephe_path(s));
}

/// C ABI: set JPL DE ephemeris file path.
#[no_mangle]
pub unsafe extern "C" fn le_set_jpl_file(fname: *const i8) {
    if fname.is_null() { return; }
    let s = unsafe { std::ffi::CStr::from_ptr(fname) }.to_str().unwrap_or("");
    with_default(|ctx| ctx.set_jpl_file(s));
}

/// C ABI: set observer position (longitude deg, latitude deg, altitude m).
#[no_mangle]
pub unsafe extern "C" fn le_set_topo(geolon: f64, geolat: f64, geoalt: f64) {
    with_default(|ctx| ctx.set_topo(geolon, geolat, geoalt));
}

/// C ABI: set sidereal (ayanamsa) mode.
#[no_mangle]
pub unsafe extern "C" fn le_set_sid_mode(sid_mode: i32, t0: f64, ayan_t0: f64) {
    with_default(|ctx| ctx.set_sid_mode(sid_mode, t0, ayan_t0));
}

/// C ABI: override Delta T (TT - UT) in seconds.
#[no_mangle]
pub unsafe extern "C" fn le_set_delta_t_user(dt: f64) {
    with_default(|ctx| ctx.set_delta_t_user(dt));
}

/// C ABI: set tidal acceleration (arcsec/cy²).
#[no_mangle]
pub unsafe extern "C" fn le_set_tid_acc(acc: f64) {
    with_default(|ctx| ctx.set_tid_acc(acc));
}

/// C ABI: load a VSOP2013 data file for a given planet (1 = Mercury .. 9 = Pluto).
#[no_mangle]
pub unsafe extern "C" fn le_set_vsop2013_file(planet: i32, path: *const i8) -> i32 {
    if path.is_null() || planet < 1 || planet > 9 { return -1; }
    let s = unsafe { std::ffi::CStr::from_ptr(path) }.to_str().unwrap_or("");
    with_default(|ctx| {
        #[cfg(feature = "analytical")]
        match ctx.set_vsop2013_file(planet as usize, s) {
            Ok(()) => 0,
            Err(_) => -1,
        }
        #[cfg(not(feature = "analytical"))]
        { let _ = s; -1 }
    })
}

/// C ABI: load ELP-MPP02 data files from a directory.
#[no_mangle]
pub unsafe extern "C" fn le_set_elpmpp02_dir(dir: *const i8) -> i32 {
    if dir.is_null() { return -1; }
    let s = unsafe { std::ffi::CStr::from_ptr(dir) }.to_str().unwrap_or("");
    with_default(|ctx| {
        #[cfg(feature = "analytical")]
        match ctx.set_elpmpp02_dir(s) {
            Ok(()) => 0,
            Err(_) => -1,
        }
        #[cfg(not(feature = "analytical"))]
        { let _ = s; -1 }
    })
}

/// C ABI: close ephemeris files and reset context.
#[no_mangle]
pub unsafe extern "C" fn le_close() {
    with_default(|ctx| {
        ctx.jpl_file_is_open = false;
        ctx.ephe_paths.clear();
        ctx.jpl_filename = None;
        #[cfg(feature = "analytical")]
        {
            ctx.vsop2013_readers.iter_mut().for_each(|r| *r = None);
            ctx.elpmpp02_reader = None;
        }
    });
}

/// C ABI: set ELP-MPP02 adjustable coefficients.
///
/// Pass a pointer to a 7-element `[f64; 7]` array:
/// `[fa, fb1, fb2, fb3, fb4, fb5, fb6]`.
/// Takes effect on the next call to le_set_elpmpp02_dir().
#[no_mangle]
pub unsafe extern "C" fn le_set_elpmpp02_paras(paras: *const f64) -> i32 {
    if paras.is_null() { return -1; }
    let arr = unsafe { std::slice::from_raw_parts(paras, 7) };
    with_default(|ctx| {
        #[cfg(feature = "analytical")]
        {
            let p = ElpParas {
                fa: arr[0], fb1: arr[1], fb2: arr[2], fb3: arr[3],
                fb4: arr[4], fb5: arr[5], fb6: arr[6],
            };
            ctx.set_elpmpp02_paras(&p);
        }
        #[cfg(not(feature = "analytical"))]
        { let _ = arr; }
    });
    0
}

/// C ABI: return library version string.
#[no_mangle]
pub unsafe extern "C" fn le_version() -> *const i8 {
    std::ffi::CStr::from_bytes_with_nul(b"0.1.0\0").unwrap().as_ptr()
}
