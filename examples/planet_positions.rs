/// Example: compute planetary positions at J2000.
/// Run with: cargo run --example planet_positions
use libre_ephemeris::calc::le_calc_ut;
use libre_ephemeris::constants;

fn main() {
    let jd_ut = 2451545.0; // J2000.0

    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL
        | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;

    let planets = [
        (constants::LE_MERCURY, "Mercury"),
        (constants::LE_VENUS, "Venus"),
        (constants::LE_MARS, "Mars"),
        (constants::LE_JUPITER, "Jupiter"),
        (constants::LE_SATURN, "Saturn"),
        (constants::LE_URANUS, "Uranus"),
        (constants::LE_NEPTUNE, "Neptune"),
    ];

    println!("Planetary positions at J2000.0 (barycentric J2000 equatorial):");
    for (ipl, name) in &planets {
        let mut xx = [0.0_f64; 24];
        let mut serr = [0_i8; 256];
        let rc = unsafe { le_calc_ut(jd_ut, *ipl, flags, xx.as_mut_ptr(), serr.as_mut_ptr()) };
        if rc == constants::LE_OK {
            let dist = (xx[0] * xx[0] + xx[1] * xx[1] + xx[2] * xx[2]).sqrt();
            println!("  {:<8}: ({:.6}, {:.6}, {:.6}) AU  dist={:.6} AU",
                     name, xx[0], xx[1], xx[2], dist);
        }
    }
}
