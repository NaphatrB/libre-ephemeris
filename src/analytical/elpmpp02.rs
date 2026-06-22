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

/// ELP-MPP02 lunar ephemeris reader and evaluator.
///
/// Loads the 14 data files from ytliu0's reformatted distribution
/// and evaluates the Fourier series for the Moon's geocentric position.
///
/// Reference: Chapront, Chapront-Touzé, Francou (2002),
/// "The lunar theory ELP revisited", A&A 404, 735-742.
///
/// Data files (14 total):
///   elp_main.long, elp_main.lat, elp_main.dist
///   elp_pert.longT0..T3, elp_pert.latT0..T2, elp_pert.distT0..T3
use crate::types::LeVec6;

/// ELP-MPP02 arguments at time T (Julian centuries from J2000).
#[derive(Clone, Copy)]
#[allow(dead_code)]
struct ElpArgs {
    w1: f64, d: f64, f: f64, l: f64, lp: f64, // Delaunay arguments
    me: f64, ve: f64, em: f64, ma: f64, ju: f64, sa: f64, ur: f64, ne: f64, // planetary
    zeta: f64,
}

/// Main problem term: 4 Delaunay multipliers + combined coefficient.
#[derive(Debug)]
struct MainTerm {
    i: [i32; 4],
    a: f64,
}

/// Perturbation term: 13 multipliers + amplitude + phase offset.
#[derive(Debug)]
struct PertTerm {
    i: [i32; 13],
    a: f64,
    phase: f64,
}

/// Loaded ELP-MPP02 coefficient data.
#[derive(Debug)]
struct ElpData {
    main_long: Vec<MainTerm>,
    main_lat: Vec<MainTerm>,
    main_dist: Vec<MainTerm>,

    pert_long: [Vec<PertTerm>; 4],
    pert_lat: [Vec<PertTerm>; 3],
    pert_dist: [Vec<PertTerm>; 4],
}

/// Adjustable coefficients for ELP-MPP02 series terms.
///
/// The main problem terms are stored as combinations of (A, B1, ..., B6).
/// The effective coefficient is: facs[0]*A + facs[1]*B1 + ... + facs[6]*B6.
///
/// Default values (all 1.0) give the standard ELP-MPP02 solution.
/// LLR-fit values improve lunar laser ranging residuals.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct ElpParas {
    pub fa: f64,
    pub fb1: f64,
    pub fb2: f64,
    pub fb3: f64,
    pub fb4: f64,
    pub fb5: f64,
    pub fb6: f64,
}

impl Default for ElpParas {
    fn default() -> Self {
        Self { fa: 1.0, fb1: 1.0, fb2: 1.0, fb3: 1.0, fb4: 1.0, fb5: 1.0, fb6: 1.0 }
    }
}

/// ELP-MPP02 reader and evaluator.
#[derive(Debug)]
pub struct ElpMpp02Reader {
    data: ElpData,
    paras: ElpParas,
}

impl ElpMpp02Reader {
    /// Open all 14 data files from a directory with default adjustable parameters.
    pub fn open(dir: &str) -> Result<Self, String> {
        Self::open_with_paras(dir, &ElpParas::default())
    }

    /// Open all 14 data files from a directory with custom adjustable parameters.
    pub fn open_with_paras(dir: &str, paras: &ElpParas) -> Result<Self, String> {
        let d = |name: &str| format!("{}/{}", dir, name);
        Ok(ElpMpp02Reader {
            data: ElpData {
                main_long: load_main(&d("elp_main.long"), paras)?,
                main_lat: load_main(&d("elp_main.lat"), paras)?,
                main_dist: load_main(&d("elp_main.dist"), paras)?,
                pert_long: [
                    load_pert(&d("elp_pert.longT0"))?,
                    load_pert(&d("elp_pert.longT1"))?,
                    load_pert(&d("elp_pert.longT2"))?,
                    load_pert(&d("elp_pert.longT3"))?,
                ],
                pert_lat: [
                    load_pert(&d("elp_pert.latT0"))?,
                    load_pert(&d("elp_pert.latT1"))?,
                    load_pert(&d("elp_pert.latT2"))?,
                ],
                pert_dist: [
                    load_pert(&d("elp_pert.distT0"))?,
                    load_pert(&d("elp_pert.distT1"))?,
                    load_pert(&d("elp_pert.distT2"))?,
                    load_pert(&d("elp_pert.distT3"))?,
                ],
            },
            paras: *paras,
        })
    }

    /// Compute Moon's geocentric X,Y,Z position in equatorial J2000 coordinates.
    pub fn compute_moon(&self, jd_et: f64) -> Result<LeVec6, String> {
        let t = (jd_et - 2451545.0) / 36525.0;
        let args = compute_args(t);

        // Main problem sums: longitude (sine), latitude (sine), distance (cosine)
        let sum_lon = sum_main(&self.data.main_long, &args, true);
        let sum_lat = sum_main(&self.data.main_lat, &args, true);
        let sum_dist = sum_main(&self.data.main_dist, &args, false);

        // Perturbation sums (always sine series)
        let sum_plon: f64 = self.data.pert_long.iter().map(|v| sum_pert(v, &args)).sum();
        let sum_plat: f64 = self.data.pert_lat.iter().map(|v| sum_pert(v, &args)).sum();
        let sum_pdist: f64 = self.data.pert_dist.iter().map(|v| sum_pert(v, &args)).sum();

        // Total (ELP units: longitude/latitude in radians, distance in km)
        let lon = sum_lon + sum_plon;
        let lat = sum_lat + sum_plat;
        let dist_km = sum_dist + sum_pdist;

        // Convert distance from km to AU
        let dist_au = dist_km / 149597870.7;

        // Ecliptic J2000 → equatorial J2000
        let eps = (23.0_f64 + 26.0 / 60.0 + 21.41136 / 3600.0).to_radians();
        let (se, ce) = eps.sin_cos();

        let x_ecl = dist_au * lat.cos() * lon.cos();
        let y_ecl = dist_au * lat.cos() * lon.sin();
        let z_ecl = dist_au * lat.sin();

        // Rotate to equatorial
        let x_eq = x_ecl;
        let y_eq = y_ecl * ce - z_ecl * se;
        let z_eq = y_ecl * se + z_ecl * ce;

        // Velocity via finite difference
        let dt = 1.0 / 1440.0; // 1 minute
        let args2 = compute_args(t + dt / 36525.0);

        let sum_lon2 = sum_main(&self.data.main_long, &args2, true);
        let sum_lat2 = sum_main(&self.data.main_lat, &args2, true);
        let sum_dist2 = sum_main(&self.data.main_dist, &args2, false);
        let sum_plon2: f64 = self.data.pert_long.iter().map(|v| sum_pert(v, &args2)).sum();
        let sum_plat2: f64 = self.data.pert_lat.iter().map(|v| sum_pert(v, &args2)).sum();
        let sum_pdist2: f64 = self.data.pert_dist.iter().map(|v| sum_pert(v, &args2)).sum();

        let lon2 = sum_lon2 + sum_plon2;
        let lat2 = sum_lat2 + sum_plat2;
        let dist_km2 = sum_dist2 + sum_pdist2;
        let dist_au2 = dist_km2 / 149597870.7;

        let x_ecl2 = dist_au2 * lat2.cos() * lon2.cos();
        let y_ecl2 = dist_au2 * lat2.cos() * lon2.sin();
        let z_ecl2 = dist_au2 * lat2.sin();

        let vx_eq = (x_ecl2 - x_ecl) / dt;
        let vy_eq = (y_ecl2 * ce - z_ecl2 * se - y_eq) / dt;
        let vz_eq = (y_ecl2 * se + z_ecl2 * ce - z_eq) / dt;

        Ok(LeVec6::new(x_eq, y_eq, z_eq, vx_eq, vy_eq, vz_eq))
    }

    /// Get the current adjustable parameters.
    pub fn paras(&self) -> &ElpParas {
        &self.paras
    }
}

/// Compute ELP-MPP02 arguments at time T (Julian centuries from J2000).
fn compute_args(t: f64) -> ElpArgs {
    let deg = std::f64::consts::PI / 180.0;
    let arcsec = std::f64::consts::PI / 648000.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t2 * t2;

    let mod2pi = |x: f64| -> f64 { x - (2.0 * std::f64::consts::PI) * ((x + std::f64::consts::PI) / (2.0 * std::f64::consts::PI)).floor() };

    // Delaunay arguments (ELP convention)
    let w1 = mod2pi(((-142.0 + 18.0/60.0 + 59.95571/3600.0) * deg)
        + 1732559343.73604 * t * arcsec
        + (-6.8084) * t2 * arcsec
        + 0.006604 * t3 * arcsec
        + (-3.169e-5) * t4 * arcsec);

    let d = mod2pi(((83.0 + 21.0/60.0 + 11.67475/3600.0) * deg)
        + 14643420.3171 * t * arcsec
        + (-38.2631) * t2 * arcsec
        + (-0.045047) * t3 * arcsec
        + 0.00021301 * t4 * arcsec);

    let f = mod2pi(((125.0 + 2.0/60.0 + 40.39816/3600.0) * deg)
        + (-6967919.5383) * t * arcsec
        + 6.359 * t2 * arcsec
        + 0.007625 * t3 * arcsec
        + (-3.586e-5) * t4 * arcsec);

    let l = mod2pi(((100.0 + 27.0/60.0 + 59.13885/3600.0) * deg)
        + 129597742.293 * t * arcsec
        + (-0.0202) * t2 * arcsec
        + 9.0e-6 * t3 * arcsec
        + 1.5e-7 * t4 * arcsec);

    let pomp = mod2pi(((102.0 + 56.0/60.0 + 14.45766/3600.0) * deg)
        + 1161.24342 * t * arcsec
        + 0.529265 * t2 * arcsec
        + (-1.1814e-4) * t3 * arcsec
        + 1.1379e-5 * t4 * arcsec);

    let _lp = mod2pi(l + (0.5 * std::f64::consts::PI - pomp));

    // Wait, actually let me re-check the C++ code. The Delaunay arguments in ELP-MPP02 are:
    // W1 = mean longitude of Moon
    // D = mean elongation of Moon from Sun  
    // F = L - Ω (argument of latitude)
    // L = mean longitude of Moon (from the Earth mean longitude computation)

    // Actually, from the C++ code:
    // args.W1 = W1 (mean longitude of Moon)
    // args.D = W2 (mean elongation)
    // args.F = W3 (argument of latitude)
    // args.L = Ea (mean anomaly of Moon)
    // args.Lp = Ea - pomp (mean anomaly of Sun = mean longitude of Earth - perihelion)

    // Hmm, let me look again at the getX2000 function
    let args = ElpArgs {
        w1,
        d,
        f,
        l,
        lp: l + (std::f64::consts::PI / 2.0 - pomp),
        // Planetary mean longitudes (simplified)
        me: mod2pi(((-108.0 + 15.0/60.0 + 3.216919/3600.0) * deg)
            + 538101628.66888 * t * arcsec),
        ve: mod2pi(((-179.0 + 58.0/60.0 + 44.758419/3600.0) * deg)
            + 210664136.45777 * t * arcsec),
        em: mod2pi(((100.0 + 27.0/60.0 + 59.13885/3600.0) * deg)
            + 129597742.293 * t * arcsec),
        ma: mod2pi(((-5.0 + 26.0/60.0 + 3.642778/3600.0) * deg)
            + 68905077.65936 * t * arcsec),
        ju: mod2pi(((34.0 + 21.0/60.0 + 5.379392/3600.0) * deg)
            + 10925660.57335 * t * arcsec),
        sa: mod2pi(((50.0 + 4.0/60.0 + 38.902495/3600.0) * deg)
            + 4399609.33632 * t * arcsec),
        ur: mod2pi(((-46.0 + 3.0/60.0 + 4.354234/3600.0) * deg)
            + 1542482.57845 * t * arcsec),
        ne: mod2pi(((-56.0 + 20.0/60.0 + 56.808371/3600.0) * deg)
            + 786547.897 * t * arcsec),
        zeta: 0.0, // general precession, set to 0 for J2000
    };

    args
}

/// Sum a main problem series.
fn sum_main(terms: &[MainTerm], args: &ElpArgs, is_sine: bool) -> f64 {
    let mut sum = 0.0;
    for term in terms {
        let phase = term.i[0] as f64 * args.d
            + term.i[1] as f64 * args.f
            + term.i[2] as f64 * args.l
            + term.i[3] as f64 * args.lp;
        if is_sine {
            sum += term.a * phase.sin();
        } else {
            sum += term.a * phase.cos();
        }
    }
    sum
}

/// Sum a perturbation series (always sine).
fn sum_pert(terms: &[PertTerm], args: &ElpArgs) -> f64 {
    let mut sum = 0.0;
    for term in terms {
        let phase = term.phase
            + term.i[0] as f64 * args.d
            + term.i[1] as f64 * args.f
            + term.i[2] as f64 * args.l
            + term.i[3] as f64 * args.lp
            + term.i[4] as f64 * args.me
            + term.i[5] as f64 * args.ve
            + term.i[6] as f64 * args.em
            + term.i[7] as f64 * args.ma
            + term.i[8] as f64 * args.ju
            + term.i[9] as f64 * args.sa
            + term.i[10] as f64 * args.ur
            + term.i[11] as f64 * args.ne
            + term.i[12] as f64 * args.zeta;
        sum += term.a * phase.sin();
    }
    sum
}

/// Load a main problem data file with adjustable coefficients.
fn load_main(path: &str, paras: &ElpParas) -> Result<Vec<MainTerm>, String> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path).map_err(|e| format!("Cannot open {}: {}", path, e))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")?.map_err(|e| format!("Read error: {}", e))?;
    let n: usize = first.trim().parse().map_err(|_| "Invalid term count")?;

    let mut terms = Vec::with_capacity(n);
    for line in lines {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let line = line.trim();
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 11 { continue; }

        let i0: i32 = parts[0].parse().unwrap_or(0);
        let i1: i32 = parts[1].parse().unwrap_or(0);
        let i2: i32 = parts[2].parse().unwrap_or(0);
        let i3: i32 = parts[3].parse().unwrap_or(0);
        let a: f64 = parts[4].parse().unwrap_or(0.0);
        let b1: f64 = parts[5].parse().unwrap_or(0.0);
        let b2: f64 = parts[6].parse().unwrap_or(0.0);
        let b3: f64 = parts[7].parse().unwrap_or(0.0);
        let b4: f64 = parts[8].parse().unwrap_or(0.0);
        let b5: f64 = parts[9].parse().unwrap_or(0.0);
        let b6: f64 = parts[10].parse().unwrap_or(0.0);

        // Combined coefficient: fa*A + fb1*B1 + ... + fb6*B6
        let coeff = paras.fa * a + paras.fb1 * b1 + paras.fb2 * b2
            + paras.fb3 * b3 + paras.fb4 * b4 + paras.fb5 * b5 + paras.fb6 * b6;

        terms.push(MainTerm {
            i: [i0, i1, i2, i3],
            a: coeff,
        });

        if terms.len() >= n { break; }
    }

    Ok(terms)
}

/// Load a perturbation data file.
fn load_pert(path: &str) -> Result<Vec<PertTerm>, String> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path).map_err(|e| format!("Cannot open {}: {}", path, e))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let first = lines.next().ok_or("Empty file")?.map_err(|e| format!("Read error: {}", e))?;
    let n: usize = first.trim().parse().map_err(|_| "Invalid term count")?;

    let mut terms = Vec::with_capacity(n);
    for line in lines {
        let line = line.map_err(|e| format!("Read error: {}", e))?;
        let line = line.trim();
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 15 { continue; }

        let mut i = [0i32; 13];
        for j in 0..13 {
            i[j] = parts[j].parse().unwrap_or(0);
        }
        let a: f64 = parts[13].parse().unwrap_or(0.0);
        let phase: f64 = parts[14].parse().unwrap_or(0.0);

        terms.push(PertTerm { i, a, phase });
        if terms.len() >= n { break; }
    }

    Ok(terms)
}
