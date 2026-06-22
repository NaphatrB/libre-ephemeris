#!/usr/bin/env python3
"""Generate comprehensive swetest oracle data for all testable features."""

import csv
import math
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
    try:
        r = subprocess.run([str(swetest)] + args, capture_output=True, text=True, timeout=30)
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
        lines.append(s)
    return lines


def sexa_to_deg(s):
    s = s.strip().replace(" ", "")
    sign = -1.0 if s.startswith("-") else 1.0
    s = s.lstrip("-")
    m = re.match(r"(\d+)°(\d+)'([\d.]+)", s)
    if not m:
        return None
    d, mi, se = float(m.group(1)), float(m.group(2)), float(m.group(3))
    return sign * (d + mi / 60.0 + se / 3600.0)


def sexa_to_rad(s):
    d = sexa_to_deg(s)
    return None if d is None else d * math.pi / 180.0


def hms_to_hours(s):
    """Convert HH:MM:SS.SS to decimal hours."""
    s = s.strip()
    m = re.match(r"(-?\d+):(\d+):([\d.]+)", s)
    if not m:
        return None
    h, mi, se = float(m.group(1)), float(m.group(2)), float(m.group(3))
    sign = -1.0 if h < 0 else 1.0
    return sign * (abs(h) + mi / 60.0 + se / 3600.0)


def main():
    swetest_path = find_swetest()
    if not swetest_path:
        print("Error: swetest binary not found.")
        sys.exit(1)

    outdir = Path(__file__).resolve().parent.parent / "tests" / "regression_data"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Using swetest: {swetest_path}")
    print(f"Output dir: {outdir}")

    # Test dates spanning analytical engine range
    jds = ["2433282.5", "2440000.5", "2445000.5", "2450000.5",
           "2451544.5", "2455000.5", "2458849.5", "2460000.5", "2469807.5"]
    planets = list(range(10))  # 0-9

    # === 1. Planet positions: equatorial J2000 cartesian (geocentric) ===
    print("1. Planet positions (equatorial J2000, geocentric)...")
    with open(outdir / "planet_positions_j2000.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 2. Planet positions: ecliptic J2000 cartesian (geocentric) ===
    print("2. Planet positions (ecliptic J2000, geocentric)...")
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

    # === 3. Planet positions: heliocentric equatorial J2000 ===
    print("3. Planet positions (heliocentric equatorial J2000)...")
    with open(outdir / "planet_positions_helio.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head", "-hel"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 4. Planet positions: barycentric equatorial J2000 ===
    print("4. Planet positions (barycentric equatorial J2000)...")
    with open(outdir / "planet_positions_bary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head", "-bary"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 5. Planet positions: equatorial of-date (no J2000) ===
    print("5. Planet positions (equatorial of-date, geocentric)...")
    with open(outdir / "planet_positions_date.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head", "-true"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 6. Planet positions: ecliptic of-date ===
    print("6. Planet positions (ecliptic of-date, geocentric)...")
    with open(outdir / "planet_positions_ecl_date.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fX", "-head", "-true"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 7. Planet positions: with all corrections (aberration, deflection, nutation) ===
    print("7. Planet positions (all corrections, equatorial J2000)...")
    with open(outdir / "planet_positions_full.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "x", "y", "z", "vx", "vy", "vz"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fx", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        w.writerow([jd, pl] + vals[:3] + [0, 0, 0])

    # === 8. Planet polar positions (longitude, latitude, distance) ===
    print("8. Planet positions (polar ecliptic J2000)...")
    with open(outdir / "planet_polar.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "lon_deg", "lat_deg", "dist_au"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fLB", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 2:
                        lon = sexa_to_deg(vals[0])
                        lat = sexa_to_deg(vals[1])
                        if lon is not None and lat is not None:
                            w.writerow([jd, pl, f"{lon:.10f}", f"{lat:.10f}", ""])

    # === 9. Planet RA/Dec (equatorial polar) ===
    print("9. Planet positions (RA/Dec equatorial J2000)...")
    with open(outdir / "planet_radec.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "ra_hours", "dec_deg", "dist_au"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-fadR", "-head"])
                if out:
                    vals = out[-1].split()
                    if len(vals) >= 3:
                        ra = float(vals[0])
                        dec = float(vals[1])
                        dist = float(vals[2])
                        w.writerow([jd, pl, f"{ra:.10f}", f"{dec:.10f}", f"{dist:.10f}"])

    # === 10. Phase, magnitude, elongation ===
    print("10. Planetary phenomena (phase, magnitude, elongation)...")
    with open(outdir / "planet_phenomena.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "phase_angle", "illuminated", "magnitude", "elongation"])
        for jd in jds:
            for pl in planets:
                out = swetest_run(swetest_path, ["-b", f"-p{pl}", f"-j{jd}", "-f+-*/", "-head"])
                if out:
                    line = out[-1]
                    # Format: "DDD°MM'SS.S" "0.xxxxx" "DDD°MM'SS.S" "DDD°MM'SS.S"
                    # Split on 2+ spaces to get individual fields
                    fields = re.split(r"\s{2,}", line)
                    if len(fields) >= 4:
                        phase_angle = sexa_to_deg(fields[0])
                        illuminated = float(fields[1].strip())
                        elongation = sexa_to_deg(fields[2])
                        magnitude = sexa_to_deg(fields[3])
                        if phase_angle is not None and elongation is not None and magnitude is not None:
                            w.writerow([jd, pl, f"{phase_angle:.6f}", f"{illuminated:.6f}",
                                       f"{magnitude:.6f}", f"{elongation:.6f}"])
                        w.writerow([jd, pl, f"{phase_angle:.6f}", f"{illuminated:.6f}",
                                   f"{magnitude:.6f}", f"{elongation:.6f}"])

    # === 11. Delta T ===
    print("11. Delta T...")
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

    # === 12. Ayanamsa ===
    print("12. Ayanamsa...")
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

    # === 13. House cusps (all systems) ===
    print("13. House cusps...")
    house_systems = ["P", "K", "E", "A", "C", "R", "M", "T", "B", "W", "X", "Y", "H", "V", "S", "N", "L", "Z"]
    locs = [("47.0", "8.5"), ("40.0", "-74.0"), ("-33.0", "151.0")]
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

    # === 14. Fixed stars ===
    print("14. Fixed stars...")
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

    # === 15. Topocentric positions ===
    print("15. Topocentric positions...")
    locs = [("47.0", "8.5", "500"), ("40.0", "-74.0", "0"), ("-33.0", "151.0", "100")]
    with open(outdir / "topocentric.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "lat", "lon", "alt", "x", "y", "z"])
        for jd in jds[:3]:
            for pl in planets[:5]:
                for lat, lon, alt in locs:
                    out = swetest_run(swetest_path, [
                        "-b", f"-p{pl}", f"-j{jd}", "-fx", "-head",
                        f"-topo{lon},{lat},{alt}"
                    ])
                    if out:
                        vals = out[-1].split()
                        if len(vals) >= 3:
                            w.writerow([jd, pl, lat, lon, alt] + vals[:3])

    # === 16. Rise/set/transit times ===
    print("16. Rise/set/transit times...")
    locs = [("47.0", "8.5"), ("40.0", "-74.0")]
    with open(outdir / "riseset.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["jd_ut", "planet", "lat", "lon", "rise_jd", "set_jd", "transit_jd"])
        for jd in jds[:3]:
            for pl in [0, 1, 2]:
                for lat, lon in locs:
                    out = swetest_run(swetest_path, [
                        "-b", f"-p{pl}", f"-j{jd}", "-rise", f"-geopos{lon},{lat}", "-head"
                    ])
                    rise = set_t = transit = None
                    for line in out:
                        if "rise" in line.lower() and "set" not in line.lower():
                            m = re.search(r"([\d.]+)", line)
                            if m: rise = m.group(1)
                        if "set" in line.lower():
                            m = re.search(r"([\d.]+)", line)
                            if m: set_t = m.group(1)
                        if "transit" in line.lower():
                            m = re.search(r"([\d.]+)", line)
                            if m: transit = m.group(1)
                    if rise or set_t or transit:
                        w.writerow([jd, pl, lat, lon, rise or "", set_t or "", transit or ""])

    print("\nDone. Generated files:")
    for p in sorted(outdir.glob("*.csv")):
        sz = p.stat().st_size
        lines = len(p.read_text().splitlines())
        print(f"  {p.name}: {sz:>6} bytes, {lines:>4} lines")


if __name__ == "__main__":
    main()
