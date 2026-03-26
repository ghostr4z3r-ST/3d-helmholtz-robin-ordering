from pathlib import Path
import sys
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/side_branches/ordering_principle/original_scripts/geometry_phase_map.py'
DEFAULT_FULL = ROOT / 'results/side_branches/ordering_principle/reference/geometry_phase_map_full.csv'
DEFAULT_SUMMARY = ROOT / 'results/side_branches/ordering_principle/reference/geometry_phase_map_summary.csv'
args = sys.argv[1:]
if '--csv' not in args: args += ['--csv', str(DEFAULT_FULL)]
if '--summary-csv' not in args: args += ['--summary-csv', str(DEFAULT_SUMMARY)]
raise SystemExit(run_historical_script(SCRIPT, args))
