#!/bin/bash
# Download ELP-MPP02 lunar ephemeris data files.
#
# Usage: ./tools/download_elpmpp02.sh [output_dir]
#   Default output: data/elpmpp02/
#
# The 14 data files are reformatted from the original ELP-MPP02 Fortran code
# by ytliu0. Source: https://github.com/ytliu0/elp-mpp02-data
#
# Files:
#   elp_main.long, elp_main.lat, elp_main.dist
#   elp_pert.longT0..T3, elp_pert.latT0..T2, elp_pert.distT0..T3

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${1:-"$SCRIPT_DIR/../data/elpmpp02"}"
mkdir -p "$OUTDIR"

BASE_URL="https://raw.githubusercontent.com/ytliu0/elp-mpp02-data/main/data"

FILES=(
    "elp_main.long"
    "elp_main.lat"
    "elp_main.dist"
    "elp_pert.longT0" "elp_pert.longT1" "elp_pert.longT2" "elp_pert.longT3"
    "elp_pert.latT0" "elp_pert.latT1" "elp_pert.latT2"
    "elp_pert.distT0" "elp_pert.distT1" "elp_pert.distT2" "elp_pert.distT3"
)

echo "Downloading ELP-MPP02 data files to $OUTDIR"
echo "Source: $BASE_URL"

for file in "${FILES[@]}"; do
    URL="${BASE_URL}/${file}"
    DEST="${OUTDIR}/${file}"
    if [ -f "$DEST" ] && [ -s "$DEST" ]; then
        echo "  $file already exists, skipping"
    else
        echo "  Downloading $file..."
        if command -v curl &>/dev/null; then
            curl -fsSL "$URL" -o "$DEST"
        elif command -v wget &>/dev/null; then
            wget -q "$URL" -O "$DEST"
        else
            echo "Error: need curl or wget"
            exit 1
        fi
        SIZE=$(stat -c%s "$DEST" 2>/dev/null || stat -f%z "$DEST" 2>/dev/null || echo "?")
        echo "    Saved: $DEST ($SIZE bytes)"
    fi
done

chmod 644 "$OUTDIR"/*
echo ""
echo "Done. Check files with: ls -lh $OUTDIR/"
ls -lh "$OUTDIR/"
