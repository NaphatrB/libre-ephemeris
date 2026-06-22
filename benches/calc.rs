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

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use libre_ephemeris::calc::le_calc_ut;
use libre_ephemeris::constants;

fn bench_planet_positions(c: &mut Criterion) {
    let jds = [
        2433282.5, 2440000.5, 2445000.5, 2450000.5, 2451544.5,
        2455000.5, 2458849.5, 2460000.5, 2469807.5,
    ];
    let planets = [2, 3, 4, 5, 6, 7, 8];
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
        | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL
        | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;

    c.bench_function("planet_positions_geocentric", |b| {
        b.iter(|| {
            for &jd in &jds {
                for &pl in &planets {
                    let mut xx = [0.0_f64; 24];
                    let mut serr = [0_i8; 256];
                    let rc = unsafe {
                        le_calc_ut(black_box(jd), black_box(pl), black_box(flags),
                                   xx.as_mut_ptr(), serr.as_mut_ptr())
                    };
                    black_box(rc);
                }
            }
        })
    });
}

fn bench_planet_positions_full(c: &mut Criterion) {
    let jds = [
        2433282.5, 2440000.5, 2445000.5, 2450000.5, 2451544.5,
        2455000.5, 2458849.5, 2460000.5, 2469807.5,
    ];
    let planets = [2, 3, 4, 5, 6, 7, 8];
    let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000;

    c.bench_function("planet_positions_full_corrections", |b| {
        b.iter(|| {
            for &jd in &jds {
                for &pl in &planets {
                    let mut xx = [0.0_f64; 24];
                    let mut serr = [0_i8; 256];
                    let rc = unsafe {
                        le_calc_ut(black_box(jd), black_box(pl), black_box(flags),
                                   xx.as_mut_ptr(), serr.as_mut_ptr())
                    };
                    black_box(rc);
                }
            }
        })
    });
}

fn bench_house_cusps(c: &mut Criterion) {
    let jds = [
        2433282.5, 2440000.5, 2445000.5, 2450000.5, 2451544.5,
        2455000.5, 2458849.5, 2460000.5, 2469807.5,
    ];
    let systems = [b'P', b'K', b'E', b'C', b'R', b'W'];

    c.bench_function("house_cusps", |b| {
        b.iter(|| {
            for &jd in &jds {
                for &hsys in &systems {
                    let mut cusps = [0.0_f64; 13];
                    let mut ascmc = [0.0_f64; 10];
                    let rc = libre_ephemeris::houses::houses(
                        black_box(jd), black_box(47.0), black_box(8.5),
                        black_box(hsys), &mut cusps, &mut ascmc,
                    );
                    black_box(rc);
                }
            }
        })
    });
}

fn bench_delta_t(c: &mut Criterion) {
    let jds = [
        2433282.5, 2440000.5, 2445000.5, 2450000.5, 2451544.5,
        2455000.5, 2458849.5, 2460000.5, 2469807.5,
    ];

    c.bench_function("delta_t", |b| {
        b.iter(|| {
            for &jd in &jds {
                let dt = unsafe { libre_ephemeris::delta_t::le_deltat(black_box(jd)) };
                black_box(dt);
            }
        })
    });
}

fn bench_coordinate_transform(c: &mut Criterion) {
    let positions = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [5.2, 0.0, 0.0],
        [-1.0, 2.0, 0.5],
    ];

    c.bench_function("equatorial_to_ecliptic", |b| {
        b.iter(|| {
            for &pos in &positions {
                let result = libre_ephemeris::transform::equatorial_to_ecliptic_cartesian(
                    black_box(pos[0]), black_box(pos[1]), black_box(pos[2]), black_box(0.4091),
                );
                black_box(result);
            }
        })
    });
}

criterion_group!(
    benches,
    bench_planet_positions,
    bench_planet_positions_full,
    bench_house_cusps,
    bench_delta_t,
    bench_coordinate_transform,
);
criterion_main!(benches);
