#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
python -m pip install -e .
./scripts/reproduce_core.sh
./scripts/main_strand/nullmodel_blindtests/reproduce_nullmodel_digest.sh
