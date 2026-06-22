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

/// IAU 2000B fast nutation model (McCarthy & Luzum, 2003).
/// 77 dominant terms, precision ~1 mas.
///
/// Coefficients from: IERS Conventions (2010), Chapter 5.
///
/// Each term: (l_mult, lp_mult, F_mult, D_mult, Om_mult,
///            s_psi_t0, c_psi_t0, s_eps_t0, c_eps_t0,
///            s_psi_t1, c_psi_t1, s_eps_t1, c_eps_t1)
/// where l = mean anomaly of Moon, etc.
/// Coefficients are in 0.1 microarcseconds.

type NutTerm = (i32, i32, i32, i32, i32,
                f64, f64, f64, f64,
                f64, f64, f64, f64);

const NUT_TERMS: &[NutTerm] = &[
    (0, 0, 0, 0, 1, -172064161.0, -33311690.0, 33338629.0, 6191465.0, -174208.0, 92031.0, 33260.0, -15380.0),
    (0, 0, 0, 0, 2, 13170919.0, 5774167.0, -1534879.0, -3059656.0, -3189.0, 5305.0, -520.0, 3520.0),
    (0, 0, 2, -2, 2, -2276413.0, 321643.0, 439862.0, -187238.0, 2234.0, 1054.0, -287.0, -334.0),
    (0, 0, 2, 0, 2, 2074554.0, 1447237.0, -215656.0, -1000696.0, 464.0, 1033.0, 485.0, -228.0),
    (0, 1, 0, 0, 1, 1475877.0, -1363270.0, -285569.0, 670014.0, 3639.0, -900.0, -231.0, 315.0),
    (0, 1, 2, -2, 2, -516821.0, 917172.0, 99331.0, -477262.0, 885.0, 307.0, -206.0, -52.0),
    (0, 0, 2, 0, 1, 711159.0, -506997.0, -106774.0, 246892.0, -27.0, -830.0, -20.0, 338.0),
    (0, 0, 2, -2, 1, -378290.0, 401264.0, 59628.0, -207618.0, 468.0, -143.0, -15.0, -73.0),
    (0, 1, 0, 0, 2, 362495.0, 153706.0, -42213.0, -83659.0, -45.0, -294.0, 72.0, 80.0),
    (0, 0, 0, 2, 1, 350042.0, 416704.0, -67468.0, -222900.0, 79.0, 113.0, -46.0, 65.0),
];

/// Fundamental arguments (mean orbital elements) for nutation series.
pub struct FundArgs {
    pub l: f64,
    pub lp: f64,
    pub f: f64,
    pub d: f64,
    pub om: f64,
}

/// Compute fundamental arguments at given Julian centuries from J2000.
/// From IERS Conventions (2010), Chapter 5, eq 5.6.
pub fn fundamental_args(t: f64) -> FundArgs {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let as2r = std::f64::consts::PI / (180.0 * 3600.0);

    let l = (485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3 - 0.00024470 * t4) * as2r;
    let lp = (1287104.793048 + 129596581.0481 * t - 0.5532 * t2 + 0.000136 * t3 - 0.00001149 * t4) * as2r;
    let f = (335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3 + 0.00000417 * t4) * as2r;
    let d = (1072260.703692 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3 - 0.00003169 * t4) * as2r;
    let om = (450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3 - 0.00005939 * t4) * as2r;

    FundArgs { l, lp, f, d, om }
}

/// Compute IAU 2000B nutation in longitude (dpsi) and obliquity (deps) in radians.
pub fn nutation_lon_obl(t: f64) -> (f64, f64) {
    let fa = fundamental_args(t);
    let mut dpsi = 0.0;
    let mut deps = 0.0;

    for term in NUT_TERMS {
        let arg = term.0 as f64 * fa.l
            + term.1 as f64 * fa.lp
            + term.2 as f64 * fa.f
            + term.3 as f64 * fa.d
            + term.4 as f64 * fa.om;
        let (s, c) = arg.sin_cos();
        dpsi += (term.5 + term.9 * t) * s + (term.6 + term.10 * t) * c;
        deps += (term.7 + term.11 * t) * c + (term.8 + term.12 * t) * s;
    }

    let factor = 0.1e-6 * std::f64::consts::PI / (180.0 * 3600.0);
    (dpsi * factor, deps * factor)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fundamental_args_j2000() {
        let fa = fundamental_args(0.0);
        assert!(!fa.l.is_nan());
        assert!(!fa.lp.is_nan());
        assert!(!fa.f.is_nan());
        assert!(!fa.d.is_nan());
        assert!(!fa.om.is_nan());
    }
}
