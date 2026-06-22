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

/// Leap second table.
///
/// Data from IERS Bulletin C (public domain).
/// Each entry: (JD of insertion, TAI-UTC offset in seconds after insertion)
/// The table covers 1972 to present.
pub const LEAP_SECONDS: &[(f64, f64)] = &[
    (2441317.5, 10.0), // 1972-01-01
    (2441499.5, 11.0), // 1972-07-01
    (2441683.5, 12.0), // 1973-01-01
    (2442048.5, 13.0), // 1974-01-01
    (2442413.5, 14.0), // 1975-01-01
    (2442778.5, 15.0), // 1976-01-01
    (2443144.5, 16.0), // 1977-01-01
    (2443509.5, 17.0), // 1978-01-01
    (2443874.5, 18.0), // 1979-01-01
    (2444239.5, 19.0), // 1980-01-01
    (2444786.5, 20.0), // 1981-07-01
    (2445151.5, 21.0), // 1982-07-01
    (2445516.5, 22.0), // 1983-07-01
    (2446247.5, 23.0), // 1985-07-01
    (2447161.5, 24.0), // 1988-01-01
    (2447892.5, 25.0), // 1990-01-01
    (2448257.5, 26.0), // 1991-01-01
    (2448804.5, 27.0), // 1992-07-01
    (2449169.5, 28.0), // 1993-07-01
    (2449534.5, 29.0), // 1994-07-01
    (2450083.5, 30.0), // 1996-01-01
    (2450630.5, 31.0), // 1997-07-01
    (2451179.5, 32.0), // 1999-01-01
    (2453736.5, 33.0), // 2006-01-01
    (2454832.5, 34.0), // 2009-01-01
    (2456109.5, 35.0), // 2012-07-01
    (2457204.5, 36.0), // 2015-07-01
    (2457754.5, 37.0), // 2017-01-01
];

/// Get the TAI-UTC offset at a given Julian day (UTC).
pub fn tai_minus_utc(jd_ut: f64) -> f64 {
    let mut offset = 10.0; // pre-1972 offset
    for &(leap_jd, ls_offset) in LEAP_SECONDS {
        if jd_ut >= leap_jd {
            offset = ls_offset;
        } else {
            break;
        }
    }
    offset
}

/// Convert UTC Julian day to TAI Julian day.
pub fn utc_to_tai(jd_utc: f64) -> f64 {
    jd_utc + tai_minus_utc(jd_utc) / 86400.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_leap_second_2020() {
        let jd = 2458849.5; // 2020-01-01
        let off = tai_minus_utc(jd);
        assert!((off - 37.0).abs() < 1e-9);
    }

    #[test]
    fn test_pre_1972() {
        let jd = 2440000.0; // 1968
        let off = tai_minus_utc(jd);
        assert_eq!(off, 10.0);
    }
}
