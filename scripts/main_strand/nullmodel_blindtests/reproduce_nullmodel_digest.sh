#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
python scripts/main_strand/nullmodel_blindtests/01_reference_digest.py
