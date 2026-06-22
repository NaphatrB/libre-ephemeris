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

pub mod reader;
pub mod interpolate;

use crate::types::LeVec6;
use reader::JplFile;

/// Compute body position using JPL DE ephemeris.
///
/// Uses the currently open JPL file from the context.
/// If no file is open, returns an error.
pub fn compute_position(jd_et: f64, ipl: i32) -> Result<LeVec6, i32> {
    crate::context::with_default(|ctx| {
        if !ctx.jpl_file_is_open || ctx.jpl_filename.is_none() {
            return Err(-1);
        }

        let mut file = JplFile::open(ctx.jpl_filename.as_ref().unwrap())?;
        let result = file.compute_body(jd_et, ipl)?;
        Ok(result)
    })
}
