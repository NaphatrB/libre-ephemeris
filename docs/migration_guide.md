# Migration Guide: Swiss Ephemeris to Libre Ephemeris

## Overview

Libre Ephemeris (`le_`) is a clean-room MIT-licensed alternative to
Swiss Ephemeris (`swe_`). The C ABI is designed for source-level
compatibility: same function argument order, same flag bit positions,
same planet numbering.

## Quick Migration

**1. Headers**

```c
// Before:
#include "swephexp.h"

// After:
#include "le_ephemeris.h"
```

**2. Prefixes**

```c
// Functions: swe_ → le_
swe_calc_ut()   → le_calc_ut()
swe_julday()    → le_julday()
swe_deltat()    → le_deltat()

// Constants: SE_ → LE_ or SEFLG_ → LE_FLG_
SEFLG_XYZ       → LE_FLG_XYZ
SEFLG_HELIO     → LE_FLG_HELIO
SE_SUN          → LE_SUN
```

**3. Ephemeris file paths**

```c
// Before: swe_set_ephe_path("path/to/ephe/files");
// After:  use de440.bsp (JPL binary) via le_set_jpl_file()
le_set_jpl_file("path/to/de440.bsp");
// Or use the analytical engine (no data files needed):
// Just call le_calc_ut() directly — VSOP87 is built-in.
```

## Complete Example

**Before (Swiss Ephemeris):**

```c
#include "swephexp.h"
#include <stdio.h>

int main() {
    char serr[256];
    double xx[6];
    int rc = swe_calc_ut(2451545.0, SE_MARS, SEFLG_XYZ, xx, serr);
    if (rc < 0) { printf("Error: %s\n", serr); return 1; }
    printf("Mars: %.6f %.6f %.6f\n", xx[0], xx[1], xx[2]);
    return 0;
}
```

**After (Libre Ephemeris):**

```c
#include "le_ephemeris.h"
#include <stdio.h>

int main() {
    char serr[256];
    double xx[24];
    int rc = le_calc_ut(2451545.0, LE_MARS, LE_FLG_XYZ, xx, serr);
    if (rc != LE_OK) { printf("Error: %s\n", serr); return 1; }
    printf("Mars: %.6f %.6f %.6f\n", xx[0], xx[1], xx[2]);
    return 0;
}
```

## Notable Differences

| Aspect | Swiss Ephemeris | Libre Ephemeris |
|--------|----------------|----------------|
| License | AGPL (commercial license required) | MIT (free for any use) |
| Data files | `.se1` (proprietary format) | JPL DE binary or none |
| Engine selection | `swe_set_ephe_path()` selects .se1 files | Flags control engine: `LE_FLG_JPLEPH`, `LE_FLG_VSOP2013`, or default VSOP87 |
| Output array | 6 doubles minimum | 24 doubles (supports multi-frame output in future) |
| Ayanamsa modes | 47 modes (same IDs) | 47 modes (same IDs) |
| House systems | Character codes ('P','K',etc.) | Character codes (same) |
| Error handling | String via `serr`, return code < 0 | String via `serr`, return code < 0 |

## Unsupported Features

These Swiss Ephemeris features are not yet implemented in Libre Ephemeris:

- SQLite-based star catalog search
- Gauquelin sector calculation
- Heliacal rise/set (phenomena)
- 1000-year Delta T extrapolation (uses known models only)
- Custom user-defined ephemeris files

## Building with Libre Ephemeris

```bash
# Build the shared library
git clone https://github.com/anomalyco/libre-ephemeris
cd libre-ephemeris
cargo build --release

# Link your program
gcc -o myapp myapp.c -I./include -L./target/release -llibre_ephemeris -lm
```

## Testing Your Migration

```bash
# Generate oracle data from your existing swetest binary
./tools/gen_regression_data.sh /path/to/swetest tests/regression_data/

# Run all regression tests
cargo test

# Validate accuracy
python3 tools/validate_accuracy.py --binary /path/to/swetest
```
