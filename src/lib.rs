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

//! Libre Ephemeris — a clean-room MIT-licensed astronomical computation library.
//!
//! Provides high-precision ephemeris calculations for planetary positions,
//! house cusps, coordinate transformations, and related astronomical data.
//!
//! ## Architecture
//!
//! The library has two ephemeris engines:
//! - **JPL**: Reads NASA DE ephemeris binary files (DE430/431/441) directly
//! - **Analytical**: Uses VSOP2013/ELP-MPP02 trigonometric series (zero-file dependency)
//!
//! The primary API is C-compatible via `extern "C"` functions with the `le_` prefix.
//! Rust consumers can also use the internal modules directly.

pub mod constants;
pub mod types;
pub mod context;
pub mod calendar;
pub mod transform;
pub mod precession;
pub mod nutation;
pub mod aberration;
pub mod deflection;
pub mod bias;
pub mod calc;
pub mod delta_t;
pub mod ayanamsa;
pub mod topocentric;
pub mod houses;
pub mod fixstar;
pub mod riseset;
pub mod data;
pub mod eclipse;
pub mod refraction;
pub mod phenomena;
pub mod heliacal;
pub mod events;

#[cfg(feature = "jpl")]
pub mod jpl;

#[cfg(feature = "analytical")]
pub mod analytical;
