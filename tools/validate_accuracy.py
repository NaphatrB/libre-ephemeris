#!/usr/bin/env python3
"""Comprehensive accuracy validation: compare le_calc vs swetest binary.

Generates N random Julian dates, runs both swetest and le_calc (via Rust test
binary or C shared library), and reports accuracy statistics for each planet.

Usage:
  # Generate comparison CSV (runs swetest only — needs manual le_calc run)
  ./tools/validate_accuracy.py --gen-only data/test_set.csv

  # Full validation with le_calc (requires compiled test binary)
  ./tools/validate_accuracy.py --binary ./swetest --output report.csv
"""

import argparse
import csv
import math
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path


# Planet names matching OE indices
PLANETS = {
    0: "Sun", 1: "Moon", 2: "Mercury", 3: "Venus", 4: "Mars",
    5: "Jupiter", 6: "Saturn", 7: "Uranus", 8: "Neptune", 9: "Pluto",
    10: "Chiron", 17: "Earth",
}

# Tolerance in arcseconds per engine type
TOLERANCE_ANALYTICAL = 5.0    # VSOP87 vs Moshier
TOLERANCE_JPL = 0.001         # JPL DE vs DE (sub-arcsecond)

# Generate test dates within analytical engine's valid range
# VSOP87 covers ±3000 years, but Moshier (swetest fallback) covers ~±50 years
# Use ±50 years from J2000 for reliable swetest comparison
DEFAULT_N_SAMPLES = 1000
ANALYTICAL_JD_MIN = 2451545.0 - 50.0 * 365.25
ANALYTICAL_JD_MAX = 2451545.0 + 50.0 * 365.25


def generate_test_set(n_samples=DEFAULT_N_SAMPLES, seed=42):
    """Generate random Julian dates and planet indices for validation."""
    rng = random.Random(seed)
    dates = set()
    while len(dates) < n_samples:
        jd = ANALYTICAL_JD_MIN + rng.random() * (ANALYTICAL_JD_MAX - ANALYTICAL_JD_MIN)
        dates.add(round(jd, 4))
    return sorted(dates)


def swetest_query(swetest_bin, jd, planet):
    """Run swetest and extract barycentric equatorial J2000 XYZ."""
    try:
        r = subprocess.run(
            [str(swetest_bin), "-b", f"-p{planet}", f"-j{jd}", "-fx", "-head"],
            capture_output=True, text=True, timeout=30,
        )
        for line in r.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith("using") or line.startswith("warning") or line.startswith("error"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                return tuple(float(p) for p in parts[:3])
    except (subprocess.TimeoutExpired, ValueError, OSError):
        pass
    return None


def angular_separation(x1, y1, z1, x2, y2, z2):
    """Compute angular separation in radians between two vectors."""
    d1 = math.sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    d2 = math.sqrt(x2 * x2 + y2 * y2 + z2 * z2)
    if d1 < 1e-15 or d2 < 1e-15:
        return 0.0
    dot = (x1 * x2 + y1 * y2 + z1 * z2) / (d1 * d2)
    return math.acos(max(-1.0, min(1.0, dot)))


def arcsec(rad):
    return rad * 206264.80624709636


def main():
    parser = argparse.ArgumentParser(description="Validate le_calc accuracy vs swetest")
    parser.add_argument("--binary", default=None, help="Path to swetest binary")
    parser.add_argument("--output", default="accuracy_report.csv", help="Output CSV path")
    parser.add_argument("--samples", type=int, default=DEFAULT_N_SAMPLES, help="Number of test dates")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--gen-only", default=None, help="Only generate test set CSV (no comparison)")
    args = parser.parse_args()

    # Generate test dates
    dates = generate_test_set(args.samples, args.seed)
    planets = sorted(PLANETS.keys())

    if args.gen_only:
        with open(args.gen_only, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["jd_ut", "planet"])
            for jd in dates:
                for pl in planets:
                    w.writerow([jd, pl])
        print(f"Generated {len(dates)} dates × {len(planets)} planets = {len(dates)*len(planets)} test cases")
        print(f"Output: {args.gen_only}")
        return

    if not args.binary:
        print("Error: --binary is required (or use --gen-only)")
        sys.exit(1)

    swetest_bin = Path(args.binary)
    if not swetest_bin.is_file():
        print(f"Error: {swetest_bin} not found")
        sys.exit(1)

    print(f"Testing {len(dates)} dates × {len(planets)} planets = {len(dates)*len(planets)} positions")
    print(f"swetest: {swetest_bin.resolve()}")

    results = []
    total = len(dates) * len(planets)
    count = 0
    failures = []

    for jd in dates:
        for pl in planets:
            count += 1
            if count % 100 == 0:
                print(f"  Progress: {count}/{total} ({100*count//total}%)")

            swetest_xyz = swetest_query(swetest_bin, jd, pl)

            if swetest_xyz is None:
                results.append({
                    "jd": jd, "planet": pl, "name": PLANETS.get(pl, "?"),
                    "sep_arcsec": None, "error": "swetest_failed",
                })
                continue

            # For now, mark result as 'swetest_only' — le_calc data
            # requires a C/Rust harness that reads this CSV and appends values.
            results.append({
                "jd": jd, "planet": pl, "name": PLANETS.get(pl, "?"),
                "sep_arcsec": None, "error": "needs_le_calc",
                "x_swetest": swetest_xyz[0],
                "y_swetest": swetest_xyz[1],
                "z_swetest": swetest_xyz[2],
            })

    # Write output
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "jd", "planet", "name", "sep_arcsec", "error",
            "x_swetest", "y_swetest", "z_swetest",
        ])
        w.writeheader()
        w.writerows(results)

    print(f"\nResults written to {args.output}")
    print(f"Total test cases: {total}")
    passed = sum(1 for r in results if r.get("sep_arcsec") is not None and r["error"] is None)
    failed = sum(1 for r in results if r.get("error") and r["error"] != "needs_le_calc")
    print(f"Passed: {passed}, Failed: {failed}, Pending le_calc: {sum(1 for r in results if r.get("error") == "needs_le_calc")}")


if __name__ == "__main__":
    main()
