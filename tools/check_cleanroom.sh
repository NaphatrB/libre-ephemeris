#!/usr/bin/env bash
# Clean-room verification script.
# Verifies that no source files from the AGPL'd Swiss Ephemeris project
# exist in this repository, and that no prohibited patterns appear in the code.
#
# Run from the project root.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"

ERRORS=0

# === Check 1: No swisseph source files exist ===
echo "[check 1] Verifying no swisseph source files..."
PROHIBITED_FILES=(
    "sweph.c"
    "sweph.h"
    "swephlib.c"
    "swephlib.h"
    "swejpl.c"
    "swejpl.h"
    "swedate.c"
    "swedate.h"
    "swehel.c"
    "swecl.c"
    "swecl.h"
    "swehouse.c"
    "swehouse.h"
    "swevents.c"
    "swevents.h"
    "swemplan.c"
    "swemmoon.c"
    "swemini.c"
    "swetest.c"
    "swephexp.h"
    "sweodef.h"
    "swemptab.h"
    "swenut2000a.h"
    "sweephe4.h"
    "swephe4.c"
    "swedll.h"
)

for f in "${PROHIBITED_FILES[@]}"; do
    if find . -name "$f" -not -path './.git/*' -not -path './tools/*' 2>/dev/null | grep -q .; then
        echo "  FAIL: Prohibited file '$f' found in repository!"
        ERRORS=$((ERRORS + 1))
    fi
done
echo "  OK"

# === Check 2: No AGPL license text ===
echo "[check 2] Verifying no AGPL license text..."
if grep -r "GNU Affero" --include="*.rs" --include="*.c" --include="*.h" --include="*.toml" --include="*.md" . 2>/dev/null | grep -v "./.git/" | grep -v "./target/" | grep -v "./tools/" | grep -q .; then
    echo "  FAIL: AGPL license text found in source files!"
    ERRORS=$((ERRORS + 1))
else
    echo "  OK"
fi

# === Check 3: No Astrodienst copyright ===
echo "[check 3] Verifying no Astrodienst copyright notices..."
if grep -ri "Astrodienst" --include="*.rs" --include="*.c" --include="*.h" --include="*.toml" --include="*.md" . 2>/dev/null \
    | grep -v "./.git/" | grep -v "./target/" | grep -v "./tools/" | grep -v "plan.md" | grep -v "No affiliation" | grep -q .; then
    echo "  FAIL: Astrodienst copyright found in source files!"
    ERRORS=$((ERRORS + 1))
else
    echo "  OK"
fi

# === Check 4: No "swe_" or "sweph_" identifiers ===
echo "[check 4] Verifying no 'swe_' or 'sweph_' identifiers in Rust source..."
if grep -rn "\bswe_\b" --include="*.rs" --include="*.h" src/ 2>/dev/null | grep -q .; then
    echo "  FAIL: 'swe_' identifier found in source code!"
    ERRORS=$((ERRORS + 1))
else
    echo "  OK"
fi

# === Check 5: All public symbols use "oe_" prefix ===
echo "[check 5] Verifying public C ABI symbols use 'oe_' prefix..."
BAD_SYMS=0
for file in $(grep -rln '#\[no_mangle\]' --include="*.rs" src/ 2>/dev/null); do
    # Read two consecutive lines, checking if #[no_mangle] is followed by a non-oe_ function
    while IFS= read -r line1 && IFS= read -r line2; do
        if echo "$line1" | grep -q '#\[no_mangle\]'; then
            if ! echo "$line2" | grep -q 'pub unsafe extern "C" fn oe_'; then
                echo "  WARNING: $file may have non-oe_ symbol"
                BAD_SYMS=$((BAD_SYMS + 1))
            fi
        fi
    done < "$file"
done
if [ "$BAD_SYMS" -eq 0 ]; then
    echo "  OK"
fi

# === Check 6: Source files have MIT license header ===
echo "[check 6] Verifying Rust source files have MIT license header..."
MISSING=0
for f in $(find src/ -name "*.rs" -not -path "./target/*"); do
    if ! head -5 "$f" | grep -qi "MIT\|Copyright"; then
        echo "  WARNING: $f missing license header"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -eq 0 ]; then
    echo "  OK"
else
    echo "  WARNING: $MISSING files without license headers"
fi

# === Summary ===
echo ""
if [ "$ERRORS" -eq 0 ]; then
    echo "CLEAN-ROOM VERIFICATION PASSED ✓"
    exit 0
else
    echo "CLEAN-ROOM VERIFICATION FAILED: $ERRORS error(s)"
    exit 1
fi
