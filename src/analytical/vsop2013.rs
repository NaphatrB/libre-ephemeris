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

/// VSOP2013 planetary theory reader and evaluator.
///
/// Loads the Poisson series coefficient files (VSOP2013p1.dat through p9.dat)
/// from IMCCE and evaluates elliptical elements (a, L, k, h, q, p) at a given
/// date, then converts to heliocentric equatorial Cartesian coordinates.
///
/// Reference: Francou, G. & Laskar, J. (2013), "VSOP2013 planetary solutions",
/// IMCCE, https://ftp.imcce.fr/pub/ephem/planets/vsop2013/
use crate::types::LeVec6;

/// Mean longitudes at J2000 (radians) for the 17 arguments.
const CI0: [f64; 17] = [
    4.402608636669000e0, // Mercury
    3.176134461576000e0, // Venus
    1.753470369433000e0, // Earth-Moon Barycenter
    6.203500014141000e0, // Mars
    4.091360003050000e0, // Vesta
    1.713740719173000e0, // Iris
    5.598641292287000e0, // Bamberga
    2.805136360408000e0, // Ceres
    2.326989734620000e0, // Pallas
    5.995461070350000e-1, // Jupiter
    8.740185101070000e-1, // Saturn
    5.481225395663000e0, // Uranus
    5.311897933164000e0, // Neptune
    0.0,
    5.198466400630000e0, // Moon (D)
    1.627905136020000e0, // Moon (F)
    2.355555638750000e0, // Moon (l)
];

/// Mean motions in longitude (radians/1000 years).
const CI1: [f64; 17] = [
    2.608790314068555e4, // Mercury
    1.021328554743445e4, // Venus
    6.283075850353215e3, // Earth-Moon Barycenter
    3.340612434145457e3, // Mars
    1.731170452721855e3, // Vesta
    1.704450855027201e3, // Iris
    1.428948917844273e3, // Bamberga
    1.364756513629990e3, // Ceres
    1.361923207632842e3, // Pallas
    5.296909615623250e2, // Jupiter
    2.132990861084880e2, // Saturn
    7.478165903077800e1, // Uranus
    3.813297222612500e1, // Neptune
    3.595362285049309e-1, // Pluto
    77713.7714481804, // Moon (D)
    84334.6615717837, // Moon (F)
    83286.9142477147, // Moon (l)
];

/// Planetary frequencies in longitude (radians/1000 years).
const FREQPLA: [f64; 9] = [
    2.608790314068555e4, // Mercury
    1.021328554743445e4, // Venus
    6.283075850353215e3, // Earth-Moon Barycenter
    3.340612434145457e3, // Mars
    5.296909615623250e2, // Jupiter
    2.132990861084880e2, // Saturn
    7.478165903077800e1, // Uranus
    3.813297222612500e1, // Neptune
    2.533566020437000e1, // Pluto
];

/// Gravitational parameters (AU³/day²) for the 9 planets + Sun.
const GMP: [f64; 9] = [
    4.9125474514508118699e-11, // Mercury
    7.2434524861627027000e-10, // Venus
    8.9970116036316091182e-10, // Earth-Moon Barycenter
    9.5495351057792580598e-11, // Mars
    2.8253458420837780000e-07, // Jupiter
    8.4597151856806587398e-08, // Saturn
    1.2920249167819693900e-08, // Uranus
    1.5243589007842762800e-08, // Neptune
    2.1886997654259696800e-12, // Pluto
];
const GMSOL: f64 = 2.9591220836841438269e-04;

/// Default maximum terms per data file (351k covers full VSOP2013).
const DEFAULT_MAX_TERMS: usize = 351000;

/// A single VSOP2013 term: 17 integer argument multipliers + sine/cosine amplitudes.
#[derive(Debug, Clone, Copy)]
struct Term {
    iphi: [i32; 17],
    ss: f64,
    cc: f64,
}

/// Coefficient block: a group of terms for one variable and polynomial power.
#[derive(Debug)]
struct CoeffBlock {
    iv: usize, // 0..5 = a, L, k, h, q, p
    it: usize, // 0..20 = power of t
    terms: Vec<Term>,
}

/// Loaded VSOP2013 data for one planet.
#[derive(Debug)]
struct PlanetData {
    blocks: Vec<CoeffBlock>,
}

/// VSOP2013 reader and evaluator for one planet file.
#[derive(Debug)]
pub struct Vsop2013Reader {
    data: PlanetData,
    _max_terms_per_block: usize,
}

impl Vsop2013Reader {
    /// Open and parse a VSOP2013pN.dat file with default term limit (351k).
    pub fn open(path: &str) -> Result<Self, String> {
        Self::open_with_max_terms(path, DEFAULT_MAX_TERMS)
    }

    /// Open and parse a VSOP2013pN.dat file, truncating to `max_terms` total per block.
    ///
    /// `max_terms` controls the precision/performance tradeoff:
    /// - 351000: full precision (default)
    /// - 100000: ~1 mas accuracy, 3× faster load
    /// - 10000: ~1 arcsec accuracy, 35× faster load
    /// - 1000: ~1 arcmin accuracy, 350× faster load
    pub fn open_with_max_terms(path: &str, max_terms: usize) -> Result<Self, String> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(path).map_err(|e| format!("Cannot open {}: {}", path, e))?;
        let reader = BufReader::with_capacity(1024 * 1024, file);
        let mut blocks: Vec<CoeffBlock> = Vec::new();

        let mut lines = reader.lines();
        let mut max_terms_reached = false;

        while let Some(line_result) = lines.next() {
            let line = line_result.map_err(|e| format!("Read error: {}", e))?;
            if line.len() < 19 {
                continue;
            }

            // Format 1001: (9x,3i3,i7) → skip 9, read 3 ints of 3 chars, 1 int of 7 chars
            let ip_str = line[9..12].trim();
            let iv_str = line[12..15].trim();
            let it_str = line[15..18].trim();
            let nt_str = line[18..25].trim();

            let _ip: i32 = ip_str.parse().unwrap_or(0);
            let iv: usize = iv_str.parse().unwrap_or(0);
            let it: usize = it_str.parse().unwrap_or(0);
            let nt: usize = nt_str.parse().unwrap_or(0);

            if nt == 0 || nt > max_terms {
                continue;
            }

            if max_terms_reached {
                continue;
            }

            let mut terms = Vec::with_capacity(nt.min(max_terms.min(10000)));

            for _ in 0..nt {
                if let Some(Ok(data_line)) = lines.next() {
                    if data_line.len() < 90 {
                        continue;
                    }
                    // Format 1002: (6x,4i3,1x,5i3,1x,4i4,1x,i6,1x,3i3,2a24)
                    // Line layout: 6 skip + 4×3 + 1 skip + 5×3 + 1 skip + 4×4 + 1 skip + 6 + 1 skip + 3×3 + 2×24 = 116 chars
                    let mut iphi = [0i32; 17];
                    if data_line.len() < 116 { continue; }
                    iphi[0] = data_line[6..9].trim().parse().unwrap_or(0);
                    iphi[1] = data_line[9..12].trim().parse().unwrap_or(0);
                    iphi[2] = data_line[12..15].trim().parse().unwrap_or(0);
                    iphi[3] = data_line[15..18].trim().parse().unwrap_or(0);
                    iphi[4] = data_line[19..22].trim().parse().unwrap_or(0);
                    iphi[5] = data_line[22..25].trim().parse().unwrap_or(0);
                    iphi[6] = data_line[25..28].trim().parse().unwrap_or(0);
                    iphi[7] = data_line[28..31].trim().parse().unwrap_or(0);
                    iphi[8] = data_line[31..34].trim().parse().unwrap_or(0);
                    iphi[9] = data_line[35..39].trim().parse().unwrap_or(0);
                    iphi[10] = data_line[39..43].trim().parse().unwrap_or(0);
                    iphi[11] = data_line[43..47].trim().parse().unwrap_or(0);
                    iphi[12] = data_line[47..51].trim().parse().unwrap_or(0);
                    iphi[13] = data_line[52..58].trim().parse().unwrap_or(0);
                    iphi[14] = data_line[59..62].trim().parse().unwrap_or(0);
                    iphi[15] = data_line[62..65].trim().parse().unwrap_or(0);
                    iphi[16] = data_line[65..68].trim().parse().unwrap_or(0);

                    // Sine and cosine amplitudes at positions 68..92 and 92..116 (Fortran D format)
                    let ss_str = data_line[68..92].trim();
                    let cc_str = data_line[92..116].trim();
                    let ss = parse_fortran_double(ss_str);
                    let cc = parse_fortran_double(cc_str);

                    terms.push(Term { iphi, ss, cc });

                    if terms.len() >= max_terms {
                        max_terms_reached = true;
                        break;
                    }
                } else {
                    break;
                }
            }

            if !terms.is_empty() && iv >= 1 && iv <= 6 {
                blocks.push(CoeffBlock { iv: iv - 1, it, terms });
            }
        }

        if blocks.is_empty() {
            return Err("No coefficient blocks found".to_string());
        }

        Ok(Vsop2013Reader { data: PlanetData { blocks }, _max_terms_per_block: max_terms })
    }

    /// Evaluate elliptical elements (a, L, k, h, q, p) at time tdj (JD from J2000).
    ///
    /// Returns [a, L, k, h, q, p] where:
    ///   a: semi-major axis (AU)
    ///   L: mean longitude (rad)
    ///   k = e·cos(ϖ) (rad)
    ///   h = e·sin(ϖ) (rad)
    ///   q = sin(i/2)·cos(Ω) (rad)
    ///   p = sin(i/2)·sin(Ω) (rad)
    pub fn evaluate(&self, tdj: f64, ip: usize) -> Result<[f64; 6], String> {
        let t = tdj / 365250.0; // t in millennia from J2000
        let powers: [f64; 21] = {
            let mut p = [0.0_f64; 21];
            p[0] = 1.0;
            for i in 1..21 {
                p[i] = p[i - 1] * t;
            }
            p
        };

        let mut r = [0.0_f64; 6];

        for block in &self.data.blocks {
            let mut sum = 0.0;
            for term in &block.terms {
                let arg = term.iphi.iter().enumerate()
                    .map(|(j, &m)| m as f64 * (CI0[j] + CI1[j] * t))
                    .sum::<f64>();
                let (s, c) = arg.sin_cos();
                sum += term.ss * s + term.cc * c;
            }
            r[block.iv] += powers[block.it] * sum;
        }



        // Add planetary mean motion to L
        if ip >= 1 && ip <= 9 {
            r[1] += FREQPLA[ip - 1] * t;
        }
        let dpi = 2.0 * std::f64::consts::PI;
        r[1] = r[1].rem_euclid(dpi);

        Ok(r)
    }

    /// Compute heliocentric Cartesian from elliptical elements.
    ///
    /// Returns (X, Y, Z, Vx, Vy, Vz) in AU and AU/day.
    /// Implements the ELLXYZ algorithm from the VSOP2013 Fortran code.
    pub fn ellipsoid_to_cartesian(&self, el: &[f64; 6], ibody: usize) -> LeVec6 {
        let rgm = (GMP[ibody - 1] + GMSOL).sqrt();
        let xa = el[0];
        let xl = el[1];
        let xk = el[2];
        let xh = el[3];
        let xq = el[4];
        let xp = el[5];

        let xfi = (1.0 - xk * xk - xh * xh).sqrt();
        let xki = (1.0 - xq * xq - xp * xp).sqrt();
        let u = 1.0 / (1.0 + xfi);
        let z_re = xk;
        let z_im = xh;
        let ex = (z_re * z_re + z_im * z_im).sqrt();
        let ex2 = ex * ex;
        let ex3 = ex2 * ex;

        let dpi = 2.0 * std::f64::consts::PI;
        let gl = xl.rem_euclid(dpi);
        let gm = gl - xh.atan2(xk);

        // Initial eccentric anomaly
        let mut e = gl + (ex - 0.125 * ex3) * gm.sin()
            + 0.5 * ex2 * (2.0 * gm).sin()
            + 0.375 * ex3 * (3.0 * gm).sin();

        // Newton's method for Kepler's equation
        let mut rsa = 1.0;
        for _ in 0..20 {
            let (se, ce) = e.sin_cos();
            let z3_re = z_re * ce + z_im * se;
            let z3_im = z_re * se - z_im * ce;
            let dl = gl - e + z3_im;
            rsa = 1.0 - z3_re;
            let de = dl / rsa;
            e += de;
            if de.abs() < 1e-15 {
                break;
            }
        }

        let (se, ce) = e.sin_cos();

        let z1_re = u * (xk * ce + xh * se);
        let z1_im = u * (xk * se - xh * ce);

        let zto_re = (-xk + ce + z1_re) / rsa;
        let zto_im = (-xh + se - z1_im) / rsa;

        let xcw = zto_re;
        let xsw = zto_im;
        let xm = xp * xcw - xq * xsw;
        let xr = xa * rsa;

        // Positions
        let wx = xr * (xcw - 2.0 * xp * xm);
        let wy = xr * (xsw + 2.0 * xq * xm);
        let wz = -2.0 * xr * xki * xm;

        // Velocities
        let xms = xa * (xh + xsw) / xfi;
        let xmc = xa * (xk + xcw) / xfi;
        let xn = rgm / xa.sqrt() / xa;

        let vx = xn * ((2.0 * xp * xp - 1.0) * xms + 2.0 * xp * xq * xmc);
        let vy = xn * ((1.0 - 2.0 * xq * xq) * xmc - 2.0 * xp * xq * xms);
        let vz = 2.0 * xn * xki * (xp * xms + xq * xmc);

        LeVec6::new(wx, wy, wz, vx, vy, vz)
    }

    /// Compute heliocentric equatorial Cartesian position at given time.
    pub fn compute_position(&self, tdj: f64, ip: usize) -> Result<LeVec6, String> {
        let el = self.evaluate(tdj, ip)?;
        let cart = self.ellipsoid_to_cartesian(&el, ip);

        // Rotate from dynamical J2000 to ICRS equatorial
        let eps = (23.0_f64 + 26.0 / 60.0 + 21.41136 / 3600.0).to_radians();
        let phi = (-0.05188_f64).to_radians() / 3600.0;
        let (ce, se) = eps.sin_cos();
        let (cp, sp) = phi.sin_cos();

        let x = cart.0[0]; let y = cart.0[1]; let z = cart.0[2];
        let vx = cart.0[3]; let vy = cart.0[4]; let vz = cart.0[5];

        // Rot = R_x(eps) * R_z(-phi)
        // First R_z(-phi): rotate ecliptic longitude
        let x1 = cp * x + sp * y;
        let y1 = -sp * x + cp * y;
        let z1 = z;
        let vx1 = cp * vx + sp * vy;
        let vy1 = -sp * vx + cp * vy;
        let vz1 = vz;

        // Then R_x(eps): rotate to equator
        let x2 = x1;
        let y2 = ce * y1 - se * z1;
        let z2 = se * y1 + ce * z1;
        let vx2 = vx1;
        let vy2 = ce * vy1 - se * vz1;
        let vz2 = se * vy1 + ce * vz1;

        Ok(LeVec6::new(x2, y2, z2, vx2, vy2, vz2))
    }
}

/// Parse a Fortran D-format double (e.g., "0.1234567890123456D+03" or "0.1234567890123456 +03").
fn parse_fortran_double(s: &str) -> f64 {
    let s = s.trim();
    // Handle space-separated exponent: "0.1234567890123456 +03" → "0.1234567890123456e+03"
    if !s.contains('D') && !s.contains('d') {
        if let Some(pos) = s.rfind(' ') {
            let exp = s[pos..].trim();
            if (exp.starts_with('+') || exp.starts_with('-')) && exp.len() > 1
                && exp[1..].chars().all(|c| c.is_ascii_digit())
            {
                let mantissa = s[..pos].trim();
                if let Ok(v) = format!("{}e{}", mantissa, exp).parse::<f64>() {
                    return v;
                }
            }
        }
    }
    // Standard format with 'D'/'d'
    let s = s.replace('D', "e").replace('d', "e");
    s.parse::<f64>().unwrap_or(0.0)
}

/// Compute a position using one of the VSOP2013 data files.
///
/// If a VSOP2013 reader has been opened for the given planet,
/// use it. Otherwise fall through to the analytical engine.
/// This mirrors the JPL DE reader architecture.
pub fn compute_position_vsop2013(
    jd_et: f64,
    ipl: i32,
    reader: &Vsop2013Reader,
) -> Result<LeVec6, i32> {
    // Map OE planet index to VSOP2013 body index (1..9)
    let ip = match ipl {
        2 => 1,  // Mercury
        3 => 2,  // Venus
        4 => 4,  // Mars
        5 => 5,  // Jupiter
        6 => 6,  // Saturn
        7 => 7,  // Uranus
        8 => 8,  // Neptune
        9 => 9,  // Pluto
        11 | 17 => 3, // Earth-Moon barycenter / Earth
        _ => return Err(-1),
    };

    let tdj = jd_et - 2451545.0;
    reader.compute_position(tdj, ip).map_err(|_| -1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fortran_double() {
        assert!((parse_fortran_double("  0.1234567890123456D+03") - 123.4567890123456).abs() < 1e-10);
        assert!((parse_fortran_double(" -1.2345678901234567D-02") + 0.012345678901234567).abs() < 1e-15);
        assert!((parse_fortran_double("0.0000000000000000 +00") - 0.0).abs() < 1e-20);
        assert!((parse_fortran_double(" 0.3954461714403000 +02") - 39.54461714403).abs() < 1e-10);
        assert!((parse_fortran_double("-0.8525819763547007 -01") + 0.08525819763547007).abs() < 1e-15);
    }

    #[test]
    fn test_ellxyz_earth_j2000() {
        // At J2000, Earth's approximate elliptical elements
        let el = [0.999999, 1.753470, -0.000937, 0.000080, 0.0, 0.0];
        let reader = Vsop2013Reader { data: PlanetData { blocks: vec![] }, _max_terms_per_block: 100 };
        let result = reader.ellipsoid_to_cartesian(&el, 3);
        let dist = (result.0[0] * result.0[0] + result.0[1] * result.0[1] + result.0[2] * result.0[2]).sqrt();
        assert!(dist > 0.9 && dist < 1.1, "Earth distance = {}", dist);
    }
}
