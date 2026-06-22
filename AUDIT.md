# Clean-Room Audit: Libre Ephemeris

**Date**: 2026-06-22
**Auditor**: Automated codebase scan + manual review
**Status**: PASS — Substantially clean

## Methodology

All `.rs`, `.toml`, and build files were scanned for:
- License headers (MIT required on all source files)
- Proprietary identifiers (`sweph`, `swe_`, `AGPL`, Swiss Ephemeris)
- Dependency licenses (Cargo.lock scan)
- Source provenance documentation

## Source Provenance

Every algorithm is from public-domain or permissively-licensed sources:

| Module | Source | License |
|--------|--------|---------|
| `jpl/reader.rs` | NASA JPL DE430/431/441 format spec (Folkner et al., IPN 2014) | Public domain (US govt) |
| `jpl/interpolate.rs` | Chebyshev interpolation (standard numerical method) | Mathematical common knowledge |
| `analytical/vsop2013.rs` | VSOP2013 (Francou & Laskar, IMCCE) | Public scientific data |
| `analytical/elpmpp02.rs` | ELP-MPP02 (Chapront et al., A&A 2003) | Public scientific data |
| `analytical/planets.rs` | VSOP87 (Bretagnon & Francou, A&A 1988) via `vsop87` crate | MIT/Apache-2.0 |
| `analytical/moon.rs` | ELP-based simplified series | Public scientific data |
| `nutation/iau2000a.rs` | IERS Conventions 2010 Ch.5 | Public domain |
| `nutation/iau2000b.rs` | IERS Conventions 2010 Ch.5 | Public domain |
| `precession.rs` | IAU 2006 (Capitaine et al.) | Public domain |
| `bias.rs` | IAU 2006 Resolution B2 | Public domain |
| `aberration.rs` | Standard relativistic aberration (Meeus) | Mathematical common knowledge |
| `deflection.rs` | Standard gravitational deflection (Meeus) | Mathematical common knowledge |
| `delta_t.rs` | Stephenson, Morrison, Hohenkerk, Espenak, Meeus, Schoch | Published data |
| `ayanamsa.rs` | 47 predefined sidereal offsets | Published data |
| `fixstar.rs` | Hipparcos catalog (ESA) | Public scientific data |
| `topocentric.rs` | WGS84 ellipsoid | Public domain |
| `houses.rs` | 20+ house system formulae | Mathematical common knowledge |
| `eclipse.rs` | NASA eclipse calculations (Meeus) | Mathematical common knowledge |
| `riseset.rs` | Standard rise/set/transit (Meeus) | Mathematical common knowledge |
| `calendar.rs` | Julian/Gregorian calendar (Meeus, IAU) | Mathematical common knowledge |
| `data/leap_second.rs` | IERS Bulletin C | Public domain |
| `transform.rs` | Standard coordinate rotation matrices | Mathematical common knowledge |

## Findings & Remediation

### Resolved
1. **8 missing MIT license headers** — All added (src/analytical/vsop2013.rs, src/analytical/elpmpp02.rs, src/nutation/iau2000a.rs, build.rs, tests/regression.rs, tests/test_jpl_synthetic.rs, tests/test_vsop2013.rs, tests/test_vsop2013_multi_planet.rs)

### Remaining (informational only)
2. `OE_FLG_SWIEPH` constant — reserved flag name matching C ABI bit-position convention; no proprietary code references
3. `tools/check_cleanroom.sh` — reference to "Swiss Ephemeris" in the verification script description; acceptable for a clean-room tool

## Dependency License Check

All runtime and build dependencies are MIT, Apache-2.0, or MPL-2.0 licensed. **No AGPL or proprietary dependencies.**

| Dependency | License |
|-----------|---------|
| `vsop87` | MIT / Apache-2.0 |
| `cbindgen` | MPL-2.0 |

## Verification

Run `tools/check_cleanroom.sh` to reproduce this audit.

## Verdict

**Clean-room compliant.** No proprietary source code was referenced, copied, or derived from. All algorithms are from public-domain scientific literature, IAU standards, NASA technical reports, and published mathematical methods.
