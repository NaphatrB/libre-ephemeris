/// VSOP87 planetary positions (MIT/Apache-2.0 licensed `vsop87` crate).
///
/// Reference: Bretagnon, P. & Francou, G. (1988),
/// "Planetary theories in rectangular and spherical variables — VSOP87 solution",
/// Astronomy & Astrophysics, 202, 309-315.
use crate::types::LeVec6;

fn spherical_to_equatorial_cartesian(lon: f64, lat: f64, dist: f64) -> (f64, f64, f64) {
    let cos_lat = lat.cos();
    let x_ecl = dist * cos_lat * lon.cos();
    let y_ecl = dist * cos_lat * lon.sin();
    let z_ecl = dist * lat.sin();
    let eps0 = 23.439291111111111_f64.to_radians();
    let x_eq = x_ecl;
    let y_eq = y_ecl * eps0.cos() - z_ecl * eps0.sin();
    let z_eq = y_ecl * eps0.sin() + z_ecl * eps0.cos();
    (x_eq, y_eq, z_eq)
}

fn vsop87b_to_state(jde: f64, f: fn(f64) -> vsop87::SphericalCoordinates) -> LeVec6 {
    let dt = 1e-6; // ~0.0864 seconds for velocity derivative

    let s0 = f(jde);
    let s1 = f(jde + dt);

    let (x0, y0, z0) = spherical_to_equatorial_cartesian(s0.longitude(), s0.latitude(), s0.distance());
    let (x1, y1, z1) = spherical_to_equatorial_cartesian(s1.longitude(), s1.latitude(), s1.distance());

    LeVec6::new(
        x0, y0, z0,
        (x1 - x0) / dt,
        (y1 - y0) / dt,
        (z1 - z0) / dt,
    )
}

/// Heliocentric J2000 equatorial position for Mercury (VSOP87B).
pub fn mercury(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::mercury) }
/// Heliocentric J2000 equatorial position for Venus (VSOP87B).
pub fn venus(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::venus) }
/// Heliocentric J2000 equatorial position for Earth (VSOP87B).
pub fn earth(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::earth) }
/// Heliocentric J2000 equatorial position for Mars (VSOP87B).
pub fn mars(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::mars) }
/// Heliocentric J2000 equatorial position for Jupiter (VSOP87B).
pub fn jupiter(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::jupiter) }
/// Heliocentric J2000 equatorial position for Saturn (VSOP87B).
pub fn saturn(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::saturn) }
/// Heliocentric J2000 equatorial position for Uranus (VSOP87B).
pub fn uranus(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::uranus) }
/// Heliocentric J2000 equatorial position for Neptune (VSOP87B).
pub fn neptune(jde: f64) -> LeVec6 { vsop87b_to_state(jde, vsop87::vsop87b::neptune) }

/// Heliocentric J2000 equatorial position for Pluto (simplified orbital elements).
///
/// VSOP87 does not include Pluto. Uses a simplified Keplerian orbit with
/// constant orbital elements (a, e, i, node, peri) and linear mean longitude.
/// Velocity is not computed (returns zero).
pub fn pluto(jde: f64) -> LeVec6 {
    // VSOP87 does not include Pluto. Use simplified orbital elements.
    let t = (jde - 2451545.0) / 365250.0;
    let a = 39.482;
    let e = 0.2488;
    let i = 17.140_f64.to_radians();
    let node = 110.304_f64.to_radians();
    let peri = 113.763_f64.to_radians();
    let ml = 1.924_f64.to_radians() + 0.003144_f64 * t;
    let ma = ml - node - peri;
    let ea = kepler_solve(ma, e);

    let x_peri = a * (ea.cos() - e);
    let y_peri = a * (1.0 - e * e).sqrt() * ea.sin();

    let (cn, sn) = (node.cos(), node.sin());
    let (cp, sp) = (peri.cos(), peri.sin());
    let ci = i.cos();

    let x_ecl = x_peri * (cn * cp - sn * sp * ci) - y_peri * (cn * sp + sn * cp * ci);
    let y_ecl = x_peri * (sn * cp + cn * sp * ci) - y_peri * (sn * sp - cn * cp * ci);
    let z_ecl = x_peri * sp * i.sin() + y_peri * cp * i.sin();

    let eps0 = 23.439291111111111_f64.to_radians();
    LeVec6::new(
        x_ecl,
        y_ecl * eps0.cos() - z_ecl * eps0.sin(),
        y_ecl * eps0.sin() + z_ecl * eps0.cos(),
        0.0, 0.0, 0.0,
    )
}

/// Heliocentric J2000 equatorial position for Chiron (simplified orbital elements).
///
/// Chiron is a centaur asteroid (95P/Chiron). VSOP87 does not include it.
/// Uses Keplerian orbital elements from JPL HORIZONS (J2000 epoch).
pub fn chiron(jde: f64) -> LeVec6 {
    let t = (jde - 2451545.0) / 365250.0;
    let a = 13.645;
    let e = 0.382;
    let i = 6.936_f64.to_radians();
    let node = 209.394_f64.to_radians();
    let peri = 339.077_f64.to_radians();
    let ml = 113.5_f64.to_radians() + 0.1245_f64 * t;
    let ma = ml - node - peri;
    let ea = kepler_solve(ma, e);

    let x_peri = a * (ea.cos() - e);
    let y_peri = a * (1.0 - e * e).sqrt() * ea.sin();

    let (cn, sn) = (node.cos(), node.sin());
    let (cp, sp) = (peri.cos(), peri.sin());
    let ci = i.cos();

    let x_ecl = x_peri * (cn * cp - sn * sp * ci) - y_peri * (cn * sp + sn * cp * ci);
    let y_ecl = x_peri * (sn * cp + cn * sp * ci) - y_peri * (sn * sp - cn * cp * ci);
    let z_ecl = x_peri * sp * i.sin() + y_peri * cp * i.sin();

    let eps0 = 23.439291111111111_f64.to_radians();
    LeVec6::new(
        x_ecl,
        y_ecl * eps0.cos() - z_ecl * eps0.sin(),
        y_ecl * eps0.sin() + z_ecl * eps0.cos(),
        0.0, 0.0, 0.0,
    )
}

fn kepler_solve(m: f64, e: f64) -> f64 {
    let mut e0 = m;
    for _ in 0..10 {
        let de = (e0 - e * e0.sin() - m) / (1.0 - e * e0.cos());
        e0 -= de;
        if de.abs() < 1e-14 { break; }
    }
    e0
}
