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

/// Synthetic JPL DE binary file tests.
///
/// Creates minimal valid JPL DE ephemeris files in a temp directory,
/// then uses the reader + interpolation pipeline to verify correctness.
use std::path::PathBuf;

const RECORD_SIZE: usize = 1656;
const RECORDS: usize = 3;

/// Write a synthetic DE header into a byte buffer.
fn write_header(data: &mut [u8], jd_start: f64, jd_end: f64) {
    let mut title = [0u8; 84];
    let msg = b"SYNTHETIC JPL DE FOR TESTING";
    title[..msg.len()].copy_from_slice(msg);
    data[..84].copy_from_slice(&title);
    data[260..268].copy_from_slice(&jd_start.to_le_bytes());
    data[268..276].copy_from_slice(&jd_end.to_le_bytes());
    data[276..280].copy_from_slice(&430i32.to_le_bytes());
}

/// Get byte offset for Mercury coefficient (set 0=X, 1=Y, 2=Z, coeff in 0..15).
fn mc(set: usize, coeff: usize) -> usize {
    // Mercury starts at byte 1656 (record 1, offset 0).
    // Each set has 15 doubles, each 8 bytes.
    RECORD_SIZE + (set * 15 + coeff) * 8
}

/// Write a Mercury coefficient into the data buffer.
fn smc(data: &mut [u8], set: usize, coeff: usize, value: f64) {
    let o = mc(set, coeff);
    data[o..o + 8].copy_from_slice(&value.to_le_bytes());
}

fn create_base() -> Vec<u8> {
    let mut data = vec![0u8; RECORDS * RECORD_SIZE];
    write_header(&mut data, 2451544.5, 2451576.5);
    data
}

fn open_synthetic(data: Vec<u8>, name: &str) -> (libre_ephemeris::jpl::reader::JplFile, PathBuf) {
    let tmp_dir = std::env::temp_dir();
    let de_path = tmp_dir.join(format!("{}.430", name));
    std::fs::write(&de_path, &data).unwrap();
    let file = libre_ephemeris::jpl::reader::JplFile::open(de_path.to_str().unwrap()).unwrap();
    (file, de_path)
}

#[test]
fn test_jpl_constant_position() {
    let mut data = create_base();
    smc(&mut data, 0, 0, 0.5);
    smc(&mut data, 1, 0, 0.3);
    smc(&mut data, 2, 0, 0.1);

    let (mut file, path) = open_synthetic(data, "jpl_const");
    let r = file.compute_body(2451545.0, 2).unwrap();

    assert!((r.0[0] - 0.5).abs() < 1e-12, "X: {}", r.0[0]);
    assert!((r.0[1] - 0.3).abs() < 1e-12, "Y: {}", r.0[1]);
    assert!((r.0[2] - 0.1).abs() < 1e-12, "Z: {}", r.0[2]);
    assert!(r.0[3].abs() < 1e-12, "Vx: {}", r.0[3]);
    assert!(r.0[4].abs() < 1e-12, "Vy: {}", r.0[4]);
    assert!(r.0[5].abs() < 1e-12, "Vz: {}", r.0[5]);
    let _ = std::fs::remove_file(&path);
}

#[test]
fn test_jpl_linear_motion() {
    let mut data = create_base();
    smc(&mut data, 0, 1, 1.0); // X = T_1(t) = t
    smc(&mut data, 1, 0, 0.5); // Y = constant 0.5

    let (mut file, path) = open_synthetic(data, "jpl_linear");

    // Segment: JD 2451544.5 to 2451576.5, 32 days.
    // t = 2*(jd - t0)/32 - 1

    // t = -1 (start): x = -1
    let r = file.compute_body(2451544.5, 2).unwrap();
    assert!((r.0[0] - (-1.0)).abs() < 1e-12, "X at t=-1: {}", r.0[0]);
    // Vx = 1.0 * (2/32) = 0.0625 AU/day
    assert!((r.0[3] - 0.0625).abs() < 1e-10, "Vx: {} != 0.0625", r.0[3]);

    // t = 0 (midpoint): x = 0
    let r = file.compute_body(2451560.5, 2).unwrap();
    assert!((r.0[0]).abs() < 1e-12, "X at t=0: {}", r.0[0]);

    // t = 1 (end): x = 1
    let r = file.compute_body(2451576.5, 2).unwrap();
    assert!((r.0[0] - 1.0).abs() < 1e-12, "X at t=1: {}", r.0[0]);

    let _ = std::fs::remove_file(&path);
}

#[test]
fn test_jpl_quadratic_motion() {
    let mut data = create_base();
    smc(&mut data, 0, 2, 1.0); // X = T_2(t) = 2t² - 1

    let (mut file, path) = open_synthetic(data, "jpl_quad");

    // t = -1: x = 2 - 1 = 1, dT_2/dt = 4t = -4 → Vx = -4*2/32 = -0.25
    let r = file.compute_body(2451544.5, 2).unwrap();
    assert!((r.0[0] - 1.0).abs() < 1e-10, "X at t=-1: {}", r.0[0]);
    assert!((r.0[3] - (-0.25)).abs() < 1e-9, "Vx at t=-1: {}", r.0[3]);

    // t = 0: x = -1, dT_2/dt = 0
    let r = file.compute_body(2451560.5, 2).unwrap();
    assert!((r.0[0] - (-1.0)).abs() < 1e-10, "X at t=0: {}", r.0[0]);
    assert!(r.0[3].abs() < 1e-12, "Vx at t=0: {}", r.0[3]);

    // t = 1: x = 1, dT_2/dt = 4 → Vx = 4*2/32 = 0.25
    let r = file.compute_body(2451576.5, 2).unwrap();
    assert!((r.0[0] - 1.0).abs() < 1e-10, "X at t=1: {}", r.0[0]);
    assert!((r.0[3] - 0.25).abs() < 1e-9, "Vx at t=1: {}", r.0[3]);

    let _ = std::fs::remove_file(&path);
}
