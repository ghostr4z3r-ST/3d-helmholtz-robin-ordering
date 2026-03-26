#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

mkdir -p results/reproduced/smoke_extended results/reproduced/side_branches

python scripts/05_center_emergence.py   --geometry cube --beta 1.0 --nx 8 --ny 8 --nz 8 --modes 10 --mode-summary 8

python scripts/06_primitive_cubic_supercell.py   --betas 1 --ncell 2 --pts-per-cell 5 --modes 10 --radius 0.35   --csv results/reproduced/smoke_extended/primitive_cubic_supercell_smoke.csv

python scripts/07_fcc_supercell_readout.py   --betas 1 --ncell 2 --pts-per-cell 5 --modes 10 --radius 0.28   --csv results/reproduced/smoke_extended/fcc_supercell_smoke.csv

python scripts/08_face_information_sidebranch.py   --betas 0,1 --nx 9 --ny 9 --nz 9 --modes 10   --csv results/reproduced/side_branches/face_information_first_test_smoke.csv
