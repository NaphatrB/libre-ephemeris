# Contributing to Libre Ephemeris

## Setup

```bash
git clone https://github.com/anomalyco/libre-ephemeris
cd libre-ephemeris
cargo build
cargo test
```

## Prerequisites

- Rust toolchain (stable, 2021 edition)
- Python 3 (for regression data generation)
- Optional: `swetest` binary (for cross-validation)

## Project Structure

```
src/
├── lib.rs              # Crate root, module declarations
├── constants.rs        # LE_* constant definitions
├── calc.rs             # Main calculation engine (le_calc)
├── context.rs          # Thread-safe state management
├── transform.rs        # Coordinate transforms
├── precession.rs       # 11 precession models
├── nutation/           # IAU 2000A + 2000B nutation
├── aberration.rs       # Annual aberration
├── deflection.rs       # Gravitational light deflection
├── bias.rs             # FK5/ICRS frame bias
├── calendar.rs         # Julian day / calendar conversions
├── delta_t.rs          # 5+ Delta T models
├── ayanamsa.rs         # 47 sidereal offset tables
├── topocentric.rs      # Observer position / parallax
├── houses.rs           # 21 house systems
├── fixstar.rs          # Hipparcos catalog reader
├── riseset.rs          # Rise/set/transit times
├── eclipse.rs          # Solar/lunar eclipse calculations
├── refraction.rs       # Atmospheric refraction
├── phenomena.rs        # Planetary magnitude, phase
├── heliacal.rs         # Heliacal phenomena
├── events.rs           # Aspect, ingress, transit search
├── analytical/         # VSOP2013 + ELP-MPP02 engines
│   ├── mod.rs
│   ├── planets.rs      # VSOP87 planet functions
│   ├── moon.rs         # Simplified ELP-MPP02
│   ├── vsop2013.rs     # VSOP2013 file reader
│   ├── elpmpp02.rs     # ELP-MPP02 file reader
│   └── series.rs       # Generic trig series evaluator
├── jpl/                # JPL DE file reader
│   ├── mod.rs
│   ├── reader.rs
│   └── interpolate.rs
└── data/
    └── leap_second.rs  # IERS leap second table
```

## Running Tests

```bash
# All tests
cargo test

# Specific test
cargo test test_regression_planet_positions_j2000

# With output
cargo test -- --nocapture

# Validation accuracy (requires swetest)
python3 tools/validate_accuracy.py --binary ./swetest --output tests/accuracy_report.csv
cargo test --test validate_accuracy -- --nocapture
```

## Generating Regression Data

```bash
python3 tools/gen_regression_data.py ./swetest tests/regression_data
```

## Code Style

- Follow existing patterns in the codebase
- No comments in code (unless absolutely necessary)
- MIT license header on all source files
- Public C ABI functions use `le_` prefix
- Rust internal API uses idiomatic types
- All functions should have doc comments for `cargo doc`

## Clean-Room Policy

This project is a clean-room implementation. All code must be derived from
public-domain sources (NASA JPL, IAU standards, VSOP2013/ELP-MPP02 published
data, Hipparcos catalog, published scientific papers). No reference to the
AGPL'd Swiss Ephemeris source code is permitted.

Run the clean-room check before submitting:
```bash
bash tools/check_cleanroom.sh
```

## Pull Request Process

1. Run `cargo test` — all tests must pass
2. Run `bash tools/check_cleanroom.sh` — must pass
3. Run `cargo fmt --check` — formatting must be clean
4. Update CHANGELOG.md if adding features or fixing bugs
5. Bump version in Cargo.toml for significant changes
