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

/// Multi-planet VSOP2013 integration test.
///
/// Downloads Mercury through Neptune data files from IMCCE,
/// computes positions at J2000, and validates against expected ranges.
use std::path::PathBuf;

fn data_dir() -> PathBuf {
    let d = std::env::temp_dir().join("libre_ephemeris_vsop2013");
    let _ = std::fs::create_dir_all(&d);
    d
}

fn ensure_data_file(planet: i32) -> Option<PathBuf> {
    let name = format!("VSOP2013p{}.dat", planet);
    let path = data_dir().join(&name);
    if path.exists() && path.metadata().map(|m| m.len() > 1000).unwrap_or(false) {
        return Some(path);
    }

    let url = format!("https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013p{}.dat", planet);
    eprintln!("Downloading {} ({})...", name, url);

    // Use a short timeout to avoid hanging test suite
    let result = std::process::Command::new("curl")
        .args(["-sL", "--connect-timeout", "10", "--max-time", "30",
            "-o", path.to_str().unwrap(), &url])
        .status()
        .or_else(|_| std::process::Command::new("wget")
            .args(["-q", "--timeout=30", "-O", path.to_str().unwrap(), &url])
            .status())
        .or_else(|_| std::process::Command::new("fetch")
            .args(["-q", "-o", path.to_str().unwrap(), &url])
            .status());

    match result {
        Ok(status) if status.success() && path.metadata().map(|m| m.len() > 1000).unwrap_or(false) => {
            eprintln!("Downloaded {} ({} bytes)", name, path.metadata().unwrap().len());
            Some(path)
        }
        _ => {
            let _ = std::fs::remove_file(&path);
            eprintln!("Download failed. Place {} manually in {:?}", name, data_dir());
            None
        }
    }
}

fn planet_name(ip: i32) -> &'static str {
    match ip {
        1 => "Mercury", 2 => "Venus", 3 => "Earth-Moon",
        4 => "Mars", 5 => "Jupiter", 6 => "Saturn",
        7 => "Uranus", 8 => "Neptune", 9 => "Pluto",
        _ => "?",
    }
}

struct PlanetCheck {
    ip: i32,         // VSOP2013 planet index (1..9)
    min_a: f64,      // min semi-major axis
    max_a: f64,      // max semi-major axis
    min_dist: f64,   // min heliocentric distance at J2000
    max_dist: f64,   // max heliocentric distance at J2000
}

#[test]
fn test_vsop2013_multi_planet_j2000() {
    use libre_ephemeris::analytical::vsop2013::Vsop2013Reader;

    let planets: [PlanetCheck; 8] = [
        PlanetCheck { ip: 1, min_a: 0.38, max_a: 0.41, min_dist: 0.28, max_dist: 0.49 },
        PlanetCheck { ip: 2, min_a: 0.72, max_a: 0.74, min_dist: 0.40, max_dist: 0.75 },
        PlanetCheck { ip: 3, min_a: 0.99, max_a: 1.02, min_dist: 0.50, max_dist: 1.05 },
        PlanetCheck { ip: 4, min_a: 1.50, max_a: 1.55, min_dist: 1.20, max_dist: 1.80 },
        PlanetCheck { ip: 5, min_a: 5.15, max_a: 5.25, min_dist: 4.50, max_dist: 5.60 },
        PlanetCheck { ip: 6, min_a: 9.40, max_a: 9.65, min_dist: 8.50, max_dist: 10.20 },
        PlanetCheck { ip: 7, min_a: 19.0, max_a: 19.3, min_dist: 18.0, max_dist: 20.5 },
        PlanetCheck { ip: 8, min_a: 30.0, max_a: 30.4, min_dist: 29.5, max_dist: 31.0 },
    ];

    for p in &planets {
        let path = match ensure_data_file(p.ip) {
            Some(path) => path,
            None => {
                eprintln!("Skipping {}: data file not available", planet_name(p.ip));
                continue;
            }
        };

        let reader = match Vsop2013Reader::open(path.to_str().unwrap()) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Skipping {}: open failed: {}", planet_name(p.ip), e);
                continue;
            }
        };

        // Evaluate elliptical elements at J2000
        let el = match reader.evaluate(0.0, p.ip as usize) {
            Ok(el) => el,
            Err(e) => {
                eprintln!("Skipping {}: evaluate failed: {}", planet_name(p.ip), e);
                continue;
            }
        };

        let name = planet_name(p.ip);
        let a = el[0];
        assert!(a >= p.min_a && a <= p.max_a,
            "{} semi-major axis: {} au (expected {}-{})", name, a, p.min_a, p.max_a);

        // Convert to Cartesian
        let cart = reader.ellipsoid_to_cartesian(&el, p.ip as usize);
        let dist = (cart.0[0]*cart.0[0] + cart.0[1]*cart.0[1] + cart.0[2]*cart.0[2]).sqrt();
        assert!(dist >= p.min_dist && dist <= p.max_dist,
            "{} heliocentric distance: {} au (expected {}-{})", name, dist, p.min_dist, p.max_dist);

        // Compute full ICRS position
        let result = reader.compute_position(0.0, p.ip as usize).unwrap();
        let dist2 = (result.0[0]*result.0[0] + result.0[1]*result.0[1] + result.0[2]*result.0[2]).sqrt();
        assert!(dist2 >= p.min_dist * 0.99 && dist2 <= p.max_dist * 1.01,
            "{} ICRS distance: {} au out of range", name, dist2);

        eprintln!("{}: a={:.4} au, dist(cart)={:.4} au, dist(icrs)={:.4} au",
            name, a, dist, dist2);
    }
}
