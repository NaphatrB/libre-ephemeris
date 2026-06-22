# Changelog

## 0.2.0 (2026-06-22)

### Added
- Refraction module (`src/refraction.rs`): Bennett (1982) formula, equatorial/horizontal coordinate conversion
- Planetary phenomena module (`src/phenomena.rs`): apparent magnitude, phase angle, illuminated fraction, elongation
- Heliacal phenomena module (`src/heliacal.rs`): arc of vision, visibility criteria, event estimation
- Event search module (`src/events.rs`): aspect finding (conjunction through biquintile), ingress detection, transit timing
- House systems: Gauquelin, Krusinski, Porphyrius (21 total systems)
- Sun geocentric regression test
- Validation accuracy test with tolerance assertions (3000" vs swetest)

### Fixed
- `mean_obliquity_iau2006`: replaced wrong `u = t/100` parameterization with correct IAU 2006 formula from Capitaine et al. (2003)
- `dt_stephenson_2016` -500 to 500 range: was using 500-1600 polynomial with wrong `u`; replaced with correct Espenak & Meeus (2006) Table 1 polynomial
- `dt_stephenson_2016` 500-1600 range: updated coefficients to match published values
- `dt_stephenson_2016` 1700-1800, 1800-1860, 1860-1900, 1900-1920, 1920-1941, 1941-1961, 1961-1986, 1986-2005 ranges: all updated with proper published coefficients
- `dt_stephenson_1984`: replaced placeholder (delegated to 1997) with actual 1984 coefficients
- `dt_espenak_meeus_2006` 1800-2000 range: split into proper sub-ranges matching the published paper
- `dt_stephenson_2016` 2000-2050 polynomial: fixed unit mismatch (was using years instead of centuries)
- Ecliptic coordinate conversion: now reads from frame-converted position instead of pre-frame-conversion position
- Frame conversion and coordinate system conversion order: frame conversion now happens before coordinate conversion
- `mean_obliquity_iau2006` return unit: was returning arcseconds instead of radians
- `houses.rs` `obliquity()`: removed double-conversion (was multiplying already-radian result by π/648000)
- `validate_accuracy.py`: fixed swetest flag order (`-b` before `-p`/`-j`), fixed date range to Moshier coverage, added error line filtering
- `check_cleanroom.sh`: allowed "No affiliation with Astrodienst" disclaimer in README.md

### Changed
- Tolerance tightened from 3600" to 3000" for regression tests
- `mean_sidereal_time_greenwich` made public in `riseset.rs`
- Version bumped to 0.2.0

## 0.1.0 (2026-06-21)

### Added
- Initial release: JPL DE reader, VSOP87 analytical engine, VSOP2013/ELP-MPP02 file readers
- All Tier 1 MVP features: planetary positions, calendar, Delta T, Ayanamsa, coordinate transforms, precession, nutation, aberration, deflection, topocentric, houses (18 systems), fixed stars
- Tier 2 features: eclipse calculations, rise/set/transit
- C ABI with `le_` prefix, cbindgen header generation
- Regression tests against swetest binary
- Clean-room verification script
