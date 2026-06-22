# Libre Ephemeris

A clean-room MIT-licensed astronomical computation library.

Compute planetary positions, Moon phases, star positions, house cusps,
eclipse data, Delta T, ayanamsa, and more. Rust implementation with a
C ABI (`le_` prefix). No affiliation with Astrodienst AG.

## License

MIT — all implementation derived from public-domain sources:
NASA JPL technical reports, IAU standards, VSOP2013/ELP-MPP02 data,
Hipparcos catalog, published scientific papers. See `AUDIT.md`.

## Quick Start

```rust
use libre_ephemeris::calc::le_calc_ut;
use libre_ephemeris::constants;

let mut xx = [0.0_f64; 24];
let mut serr = [0_i8; 256];
let jd = 2451545.0; // J2000.0
let flags = constants::LE_FLG_XYZ | constants::LE_FLG_J2000
    | constants::LE_FLG_NOABERR | constants::LE_FLG_NOGDEFL
    | constants::LE_FLG_NOBIRR | constants::LE_FLG_NONUT;

let rc = unsafe {
    le_calc_ut(jd, constants::LE_MARS, flags,
               xx.as_mut_ptr(), serr.as_mut_ptr())
};
assert_eq!(rc, 0);
println!("Mars: ({:.6}, {:.6}, {:.6}) AU", xx[0], xx[1], xx[2]);
```

## API Reference

### Ephemeris Calculation

| Function | Purpose |
|----------|---------|
| `le_calc(tjd, ipl, iflag, xx, serr)` | Main calculation (ET input) |
| `le_calc_ut(tjd_ut, ipl, iflag, xx, serr)` | Main calculation (UT input) |
| `le_sol_eclipse_how(tjd, iflag, attr, serr)` | Solar eclipse data |
| `le_lun_eclipse_how(tjd, iflag, attr, serr)` | Lunar eclipse data |

### Calendar

| Function | Purpose |
|----------|---------|
| `le_julday(year, month, day, gregflag)` | Date → Julian day number |
| `le_revjul(jd, gregflag, year, month, day)` | Julian day → date |
| `le_day_of_week(jd)` | Day of week (0=Mon..6=Sun) |

### Coordinate Transformation

| Function | Purpose |
|----------|---------|
| `le_cotrans(x, y, z, eps, xout, yout, zout)` | Equatorial ↔ ecliptic |
| `le_cotrans_sp(x, y, z, eps, xout, yout, zout)` | With speed |
| `le_split_deg(ddeg, roundflag, deg, min, sec, dsecfrac)` | Degrees → DMS |
| `le_csnorm(d)` | Normalize to [0, 360) |
| `le_difdeg2n(d1, d2)` | Minimum angle difference |

### Delta T

| Function | Purpose |
|----------|---------|
| `le_deltat(tjd)` | Delta T (TT − UT) in seconds |
| `le_deltat_ex(tjd, iflag, serr)` | With extended interface |

### Ayanamsa

| Function | Purpose |
|----------|---------|
| `le_get_ayanamsa(tjd_et)` | Compute ayanamsa |
| `le_get_ayanamsa_ex(tjd_et, iflag, serr)` | With extended interface |
| `le_get_ayanamsa_name(sidmode)` | Get ayanamsa name by ID |

### Houses

| Function | Purpose |
|----------|---------|
| `le_houses(tjd_ut, lat, lon, hsys, cusps, ascmc)` | Compute houses |
| `le_houses_ex(tjd_ut, iflag, lat, lon, hsys, cusps, ascmc, serr)` | With error |
| `le_house_pos(armc, geolat, eps, hsys, xpin, serr)` | House position of ecliptic point |
| `le_house_name(hsys)` | House system name |

### Fixed Stars

| Function | Purpose |
|----------|---------|
| `le_fixstar(name, jd_et, iflag, xx, serr)` | Star position at date |
| `le_star_count()` | Number of stars (200) |
| `le_star_data(index, star)` | Star data by index |

### Rise/Set/Transit

| Function | Purpose |
|----------|---------|
| `le_rise_trans(tjd_ut, ipl, iflag, rsmi, lon, lat, attr, serr)` | Rise/set/transit times |

### Context Configuration

| Function | Purpose |
|----------|---------|
| `le_set_ephe_path(path)` | Set ephemeris search paths |
| `le_set_jpl_file(fname)` | Set JPL DE file path |
| `le_set_topo(lon, lat, alt)` | Set observer position |
| `le_set_sid_mode(mode, t0, ayan_t0)` | Set ayanamsa mode |
| `le_set_delta_t_user(dt)` | Override Delta T |
| `le_set_tid_acc(acc)` | Set tidal acceleration |
| `le_set_vsop2013_file(planet, path)` | Load VSOP2013 data |
| `le_set_elpmpp02_dir(dir)` | Load ELP-MPP02 data |
| `le_set_elpmpp02_paras(paras)` | Set ELP adjustable params |
| `le_close()` | Reset context |
| `le_version()` | Library version string |

### Planet Index Constants

| Constant | Value | Body |
|----------|-------|------|
| `LE_SUN` | 0 | Sun |
| `LE_MOON` | 1 | Moon |
| `LE_MERCURY` | 2 | Mercury |
| `LE_VENUS` | 3 | Venus |
| `LE_MARS` | 4 | Mars |
| `LE_JUPITER` | 5 | Jupiter |
| `LE_SATURN` | 6 | Saturn |
| `LE_URANUS` | 7 | Uranus |
| `LE_NEPTUNE` | 8 | Neptune |
| `LE_PLUTO` | 9 | Pluto |
| `LE_CHIRON` | 10 | Chiron |
| `LE_MEAN_BARY` | 11 | Earth-Moon barycenter |
| `LE_EARTH` | 17 | Earth |

### Flag Constants

| Flag | Value | Meaning |
|------|-------|---------|
| `LE_FLG_VSOP2013` | 0x0001 | Use VSOP2013 engine |
| `LE_FLG_JPLEPH` | 0x0002 | Use JPL DE engine |
| `LE_FLG_HELIO` | 0x0010 | Heliocentric output |
| `LE_FLG_BARYHEL` | 0x0020 | Barycentric output |
| `LE_FLG_TOPOCTR` | 0x0040 | Topocentric output |
| `LE_FLG_XYZ` | 0x0080 | Cartesian output |
| `LE_FLG_SPEED` | 0x0100 | Include velocity |
| `LE_FLG_NOABERR` | 0x0200 | Skip aberration |
| `LE_FLG_NOGDEFL` | 0x0400 | Skip deflection |
| `LE_FLG_NOBIRR` | 0x0800 | Skip frame bias |
| `LE_FLG_J2000` | 0x1000 | J2000 frame output |
| `LE_FLG_SIDEREAL` | 0x4000 | Sidereal output |
| `LE_FLG_ICRS` | 0x8000 | ICRS frame output |
| `LE_FLG_EQUATORIAL` | 0x10000 | Equatorial output |
| `LE_FLG_ECLIPTIC` | 0x20000 | Ecliptic output |
| `LE_FLG_NONUT` | 0x100000 | Skip nutation |

### Error Codes

| Code | Value | Meaning |
|------|-------|---------|
| `LE_OK` | 0 | Success |
| `LE_ERR` | −1 | General error |
| `LE_ERR_INVALID_PARAMS` | −2 | Invalid parameters |
| `ERR_FILE_NOT_FOUND` | −3 | Data file not found |
| `ERR_OUT_OF_RANGE` | −4 | Parameter out of range |
| `ERR_NO_EPHEMERIS` | −5 | No ephemeris loaded |
| `ERR_IO` | −6 | I/O error |
| `ERR_ENGINE` | −7 | Engine computation failed |
| `ERR_INVALID_PLANET` | −8 | Invalid planet index |
| `ERR_INVALID_FLAG` | −9 | Invalid flag combination |
| `ERR_NOT_IMPLEMENTED` | −11 | Feature not implemented |

## C Example

```c
#include "le_ephemeris.h"
#include <stdio.h>

int main() {
    double xx[24];
    char serr[256];
    int rc = le_calc_ut(2451545.0, LE_MARS,
        LE_FLG_XYZ | LE_FLG_J2000 | LE_FLG_NOABERR
        | LE_FLG_NOGDEFL | LE_FLG_NOBIRR | LE_FLG_NONUT,
        xx, serr);
    if (rc == LE_OK) {
        printf("Mars: %.6f %.6f %.6f AU\n", xx[LE_X], xx[LE_Y], xx[LE_Z]);
    }
    return 0;
}
```

Compile with:
```bash
gcc -o example example.c -llibre_ephemeris -lm
```

## Building

```bash
cargo build --release        # Static + shared library
cargo build --features jpl   # With JPL DE support
```

## Testing

```bash
cargo test                   # Unit + integration tests
cargo test -- --ignored      # Data-dependent tests (VSOP2013 download)
./tools/gen_regression_data.sh ./swetest tests/regression_data/
cargo test                   # Now includes swetest cross-validation
```

## Migration from Swiss Ephemeris

The `le_` C ABI is designed to be a drop-in replacement at the
calling convention level. Most flag constants, planet indices,
and function names use the same numeric values.

**Key differences:**
- No `swe_` prefix — uses `le_` prefix
- No proprietary `.se1` data files — uses JPL DE binary format
- Analytical engine (VSOP2013 + ELP-MPP02) needs zero data files
- Precession, nutation, and frame bias models are independently selectable
- MIT licensed — no AGPL restrictions on commercial use

**Migration steps:**
1. Replace `#include "swephexp.h"` → `#include "le_ephemeris.h"`
2. Replace `swe_` → `le_` prefix on all function calls
3. Replace `SE_` → `LE_` prefix on all constants
4. Use `le_calc_ut()` instead of `swe_calc_ut()`
5. Use `le_set_jpl_file()` instead of `swe_set_ephe_path()` for DE files
6. Analytical engine works without data files — just call `le_calc_ut()`

## Ephemeris Engines

| Engine | Precision | Data | Use Case |
|--------|-----------|------|----------|
| VSOP87 (default) | ~1" | None | Fast, always available |
| VSOP2013 | ~0.1" | 500 MB download | High-precision analytical |
| ELP-MPP02 (Moon) | ~0.01" | 2 MB download | Accurate Moon |
| JPL DE430/431/441 | ~0.001" | 2 GB file | Research-grade |

## Data Files

```bash
# Download VSOP2013 for high-precision planets (~500 MB total)
./tools/download_vsop2013.sh data/vsop2013/

# Download ELP-MPP02 for accurate Moon (~2 MB)
./tools/download_elpmpp02.sh data/elpmpp02/

# Generate regression oracle from swetest binary
./tools/gen_regression_data.sh ./swetest tests/regression_data/
```

## Project Status

- [x] Planetary positions (VSOP87, VSOP2013, JPL DE)
- [x] Lunar positions (simplified, ELP-MPP02)
- [x] Coordinate transformations
- [x] Precession (11 models)
- [x] Nutation (IAU 2000A, 2000B, 1980, Woolard)
- [x] Aberration, light deflection
- [x] Frame bias (IAU 2000, 2006)
- [x] Delta T (5 models)
- [x] Ayanamsa (47 modes)
- [x] House cusps (20+ systems)
- [x] Fixed stars (200-star catalog)
- [x] Rise/set/transit times
- [x] Solar/lunar eclipses
- [x] C ABI (35 exported functions)
- [x] Clean-room audit (PASS)
