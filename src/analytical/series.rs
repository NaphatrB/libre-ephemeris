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

/// Generic trigonometric series evaluator for VSOP2013 and ELP-MPP02.
///
/// Both VSOP2013 (planets) and ELP-MPP02 (Moon) provide positions as
/// sums of periodic terms:
///   Σ A * cos(B + C * t)
///
/// where A, B, C are the amplitude, phase, and frequency coefficients
/// from the published series, and t is time in millennia from J2000.

/// A single periodic term in a trig series.
#[derive(Debug, Clone, Copy)]
pub struct SeriesTerm {
    pub amp: f64,     // amplitude A
    pub phase: f64,   // phase B (in radians)
    pub freq: f64,    // frequency C (in rad/millennium)
}

/// Evaluate a trigonometric series at time t (in millennia from J2000).
///
/// Returns Σ A * cos(B + C * t)
pub fn evaluate_series(terms: &[SeriesTerm], t: f64) -> f64 {
    let mut sum = 0.0;
    for term in terms {
        sum += term.amp * (term.phase + term.freq * t).cos();
    }
    sum
}

/// Evaluate a trigonometric series and its derivative at time t.
///
/// Returns (value, derivative w.r.t. t).
pub fn evaluate_series_with_deriv(terms: &[SeriesTerm], t: f64) -> (f64, f64) {
    let mut val = 0.0;
    let mut deriv = 0.0;
    for term in terms {
        let arg = term.phase + term.freq * t;
        let (s, c) = arg.sin_cos();
        val += term.amp * c;
        deriv -= term.amp * term.freq * s;
    }
    (val, deriv)
}

/// Evaluate VSOP2013 style: elliptical coordinates (L, B, R) then transform.
/// L = longitude, B = latitude, R = radius.
pub fn evaluate_vsop(
    terms_l: &[SeriesTerm],
    terms_b: &[SeriesTerm],
    terms_r: &[SeriesTerm],
    t: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    // Evaluate series for each coordinate
    let (lon, lon_dot) = evaluate_series_with_deriv(terms_l, t);
    let (lat, lat_dot) = evaluate_series_with_deriv(terms_b, t);
    let (rad, rad_dot) = evaluate_series_with_deriv(terms_r, t);

    (lon, lat, rad, lon_dot, lat_dot, rad_dot)
}

/// Precomputed series terms for a planet.
/// In a full implementation, these would be generated from the VSOP2013/ELP
/// data files at build time. Here we provide representative terms for
/// verification purposes. The full series has thousands of terms.
pub struct PlanetSeries {
    pub name: &'static str,
    pub l: &'static [SeriesTerm],
    pub b: &'static [SeriesTerm],
    pub r: &'static [SeriesTerm],
}
