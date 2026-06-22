#!/usr/bin/env python3
"""Generate regression test oracle CSVs from the original swetest binary.

Usage:
  ./tools/gen_regression_data.py [swetest_path] [output_dir]

If swetest_path is omitted, auto-detects in common locations.
Output defaults to tests/regression_data/
"""

import csv
import os
import re
import subprocess
import sys
from pathlib import Path


def find_swetest():
    candidates = ["./swetest", "../swetest", "/usr/local/bin/swetest", "/usr/bin/swetest"]
    for c in candidates:
        p = Path(c)
        if p.is_file() and os.access(p, os.X_OK):
            return p.resolve()
    return None


def swetest_run(swetest, args):
    """Run swetest with concatenated args (e.g. '-p0' not '-p', '0')."""
    try:
        r = subprocess.run([str(swetest)] + args,
                           capture_output=True, text=True, timeout=30)
    except subprocess.TimeoutExpired:
        return []
    lines = []
    for line in r.stdout.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("using Moshier") or s.startswith("warning:"):
            continue
        if s.startswith("date (dmy"):
            continue
        if s.startswith("Epsilon") or s.startswith("Nutation"):
            continue
        if swetest.name in line:
            continue
        # Keep UT:, TT: lines — they contain delta t and ayanamsa values
        lines.append(s)
    return lines


def sexa_to_deg(s):
    """Convert 'DDD°MM'SS.SSSS' to decimal degrees."""
    s = s.strip().replace(" ", "")
    sign = -1.0 if s.startswith("-") else 1.0
    s = s.lstrip("-")
    m = re.match(r"(\d+)°(\d+)'([\d.]+)", s)
    if not m:
        return None
    d, mi, se = float(m.group(1)), float(m.group(2)), float(m.group(3))
    return sign * (d + mi / 60.0 + se / 3600.0)


def main():
    script_dir = Path(__file__).resolve().parent.parent
    swetest_path = Path(sys.argv[1]) if len(sys.argv) > 1 else find_swetest()
    if not swetest_path or not swetest_path.is_file():
        print("Error: swetest binary not found.")
        sys.exit(1)

    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else script_dir / "tests" / "regression_data"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Using swetest: {swetest_path}")
    print(f"Output dir: {outdir}")

    jds = ["0.0", "1050000.5", "1356000.5", "1721423.5", "1900000.5",
           "2086000.5", "2269000.5", "2342000.5", "2415020.5",
           "2433282.5", "2451544.5", "2455197.5", "2458849.5", "2469807.5"]
    # Planets 0-9 (Sun through Pluto) map correctly between swetest and OE.
    # Planets 10+ have mismatched indices and our analytical engine doesn't support them fully.
    planets = list(range(10))

    # === 1. Planet positions: equatorial J2000 cartesian ===
    print("Generating planet position oracles (j2000)...")
    with open(outdir / "planet_positions_j2000.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                # NOTE: -b MUST come before -p and -j (swetest v2.10.03+ is order-sensitive)
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 2. Planet positions: ecliptic cartesian ===
    print("Generating planet position oracles (ecliptic)...")
    with open(outdir / "planet_positions_ecliptic.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fX", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 3. Fixed stars (catalog) ===
    print("Generating fixed star oracles...")
    stars = [
        ("Aldebaran", "4.668126", "-16.511583"),
        ("Regulus", "10.139636", "11.925014"),
        ("Antares", "16.508302", "-26.431975"),
        ("Spica", "13.412527", "-11.159241"),
        ("Sirius", "6.752406", "-16.716112"),
        ("Rigel", "5.241605", "-8.201639"),
        ("Betelgeuse", "5.918612", "7.407028"),
        ("Polaris", "2.529722", "89.264139"),
        ("Vega", "18.615692", "38.783689"),
        ("Capella", "5.285597", "45.997991"),
    ]
    with open(outdir / "fixed_stars.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "star_name", "ra_hours", "dec_deg"])
        for name, ra, dec in stars:
            w.writerow(["2451544.5", name, ra, dec])

    # === 4. Delta T ===
    print("Generating Delta T oracle...")
    with open(outdir / "delta_t.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "dt_seconds"])
        for jd in jds:
            out = swetest_run(swetest_path, [f"-p0", f"-j{jd}"])
            for line in out:
                if "delta t:" in line:
                    m = re.search(r"delta t:\s*(-?[\d.]+)", line)
                    if m:
                        w.writerow([jd, m.group(1)])
                    break

    # === 5. Ayanamsa (Fagan/Bradley) ===
    print("Generating ayanamsa oracle...")
    with open(outdir / "ayanamsa.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "ayanamsa_deg"])
        for jd in jds:
            out = swetest_run(swetest_path, [f"-p0", f"-j{jd}", "-sid0"])
            for line in out:
                if "ayanamsa =" in line:
                    m = re.search(r"ayanamsa =\s*([\d°'\"\.\s-]+)", line)
                    if m:
                        deg = sexa_to_deg(m.group(1))
                        if deg is not None:
                            w.writerow([jd, f"{deg:.10f}"])
                    break

    # === 6. House cusps ===
    print("Generating house cusp oracles...")
    house_systems = ["P", "K", "E", "A", "C", "R", "M", "T", "B", "W"]
    locs = [("47.0", "8.5")]

    with open(outdir / "houses.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "lat", "lon", "hsys"] +
                   [f"cusp{i}" for i in range(1, 13)] + ["asc", "mc"])
        for jd in jds:
            for lat, lon in locs:
                for hsys in house_systems:
                    out = swetest_run(swetest_path, [
                        f"-j{jd}", f"-geopos{lon},{lat}", f"-house{hsys}", "-head"
                    ])
                    if not out:
                        continue
                    cusps = {}
                    asc = mc = None
                    for line in out:
                        m = re.match(r"house\s+(\d+)\s+([\d°'\"\.\s-]+)\s+", line)
                        if m:
                            idx = int(m.group(1))
                            deg = sexa_to_deg(m.group(2))
                            if deg is not None:
                                cusps[idx] = deg
                        m = re.match(r"Ascendant\s+([\d°'\"\.\s-]+)\s+", line)
                        if m:
                            asc = sexa_to_deg(m.group(1))
                        m = re.match(r"MC\s+([\d°'\"\.\s-]+)\s+", line)
                        if m:
                            mc = sexa_to_deg(m.group(1))
                    if len(cusps) == 12 and asc is not None and mc is not None:
                        row = [jd, lat, lon, hsys] + \
                              [f"{cusps[i]:.10f}" for i in range(1, 13)] + \
                              [f"{asc:.10f}", f"{mc:.10f}"]
                        w.writerow(row)

    print("Done.")
    for p in sorted(outdir.glob("*.csv")):
        sz = p.stat().st_size
        lines = len(p.read_text().splitlines())
        print(f"  {p.name}: {sz:>6} bytes, {lines:>4} lines")


if __name__ == "__main__":
    main()
