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

/// Integration test for VSOP2013 reader using Pluto (smallest data file, ~13 MB).
///
/// Downloads the Pluto data file from IMCCE if not cached locally,
/// then validates a computed position.
use std::path::PathBuf;
use std::process::Command;

fn data_dir() -> PathBuf {
    let d = std::env::temp_dir().join("libre_ephemeris_vsop2013");
    let _ = std::fs::create_dir_all(&d);
    d
}

fn ensure_data_file() -> Option<PathBuf> {
    let path = data_dir().join("VSOP2013p9.dat");
    if path.exists() && path.metadata().map(|m| m.len() > 1000).unwrap_or(false) {
        return Some(path);
    }

    let url = "https://ftp.imcce.fr/pub/ephem/planets/vsop2013/solution/VSOP2013p9.dat";
    eprintln!("Downloading VSOP2013p9.dat (Pluto, ~13 MB) from {url}...");

    let result = Command::new("curl")
        .args(["-sL", "--connect-timeout", "10", "--max-time", "30",
            "-o", path.to_str().unwrap(), url])
        .status()
        .or_else(|_| Command::new("wget")
            .args(["-q", "--timeout=30", "-O", path.to_str().unwrap(), url])
            .status())
        .or_else(|_| Command::new("fetch")
            .args(["-q", "-o", path.to_str().unwrap(), url])
            .status());

    match result {
        Ok(status) if status.success() && path.metadata().map(|m| m.len() > 1000).unwrap_or(false) => {
            eprintln!("Downloaded ({})", path.metadata().unwrap().len());
            Some(path)
        }
        _ => {
            let _ = std::fs::remove_file(&path);
            eprintln!("Download failed. Place VSOP2013p9.dat manually in {:?}", data_dir());
            None
        }
    }
}

#[test]
fn test_vsop2013_pluto_j2000() {
    let path = match ensure_data_file() {
        Some(p) => p,
        None => {
            eprintln!("Skipping: VSOP2013p9.dat not available.");
            return;
        }
    };

    let reader = libre_ephemeris::analytical::vsop2013::Vsop2013Reader::open(
        path.to_str().unwrap()
    ).expect("Failed to open VSOP2013p9.dat");

    let el = reader.evaluate(0.0, 9).expect("Failed to evaluate Pluto");

    // Pluto semi-major axis ~39.5 AU
    assert!(el[0] > 35.0 && el[0] < 45.0, "Pluto a = {} au", el[0]);

    // Convert to Cartesian
    let cart = reader.ellipsoid_to_cartesian(&el, 9);
    let d = (cart.0[0]*cart.0[0] + cart.0[1]*cart.0[1] + cart.0[2]*cart.0[2]).sqrt();
    assert!(d > 25.0 && d < 55.0, "Pluto dist = {} au", d);

    // Full ICRS position
    let result = reader.compute_position(0.0, 9).unwrap();
    let d2 = (result.0[0]*result.0[0] + result.0[1]*result.0[1] + result.0[2]*result.0[2]).sqrt();
    assert!(d2 > 25.0 && d2 < 55.0, "Pluto ICRS dist = {} au", d2);

    eprintln!("Pluto at J2000: ({:.6}, {:.6}, {:.6}) au",
        result.0[0], result.0[1], result.0[2]);
}
