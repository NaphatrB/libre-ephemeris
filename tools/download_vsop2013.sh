#!/bin/bash
# Download VSOP2013 planetary theory data files from IMCCE.
#
# Usage: ./tools/download_vsop2013.sh [output_dir]
#   Default output: data/vsop2013/
#
# Downloads 9 files: VSOP2013p1.dat through VSOP2013p9.dat
# Total size: ~500 MB

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="${1:-"$SCRIPT_DIR/../data/vsop2013"}"
mkdir -p "$OUTDIR"

BASE_URL="https://ftp.imcce.fr/pub/ephem/planets/vsop2013"

echo "Downloading VSOP2013 data files to $OUTDIR"
echo "Source: $BASE_URL"

for i in $(seq 1 9); do
    FILE="VSOP2013p${i}.dat"
    URL="${BASE_URL}/${FILE}"
    DEST="${OUTDIR}/${FILE}"
    if [ -f "$DEST" ] && [ -s "$DEST" ]; then
        echo "  $FILE already exists, skipping"
    else
        echo "  Downloading $FILE..."
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

echo ""
echo "Done. Check files with: ls -lh $OUTDIR/"
ls -lh "$OUTDIR/"
