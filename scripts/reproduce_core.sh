#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="$ROOT/results/reproduced"
mkdir -p "$OUT"

python scripts/01_cube_solver_basics.py --n 10 --beta 1.0 --modes 8 --csv-dir "$OUT/cube"
python scripts/02_geometry_center_readout.py --beta 1.0 --csv "$OUT/geometry_center_beta1.csv"
python scripts/03_tetrahedron_control.py --beta 1.0 --csv "$OUT/tetrahedron_beta1.csv"
python scripts/04_hexahedron_shear_study.py --betas 1,4.8 --csv "$OUT/hexa_shear_beta1_4p8.csv"

echo "Reproduced core outputs in $OUT"
