#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
python -m pip install -e .
python scripts/side_branches/face_diagnosis/01_face_information_first_test.py
