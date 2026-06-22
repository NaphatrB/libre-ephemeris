#!/usr/bin/env bash
# Generate regression test oracle data from the original compiled swetest binary.
#
# Usage:
#   ./tools/gen_regression_data.sh /path/to/swetest [output_dir]
#
# Delegates to the Python implementation for robust CSV generation.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/gen_regression_data.py" "$@"
