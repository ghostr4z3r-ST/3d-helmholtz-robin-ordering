from pathlib import Path
import sys
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/side_branches/ordering_principle/original_scripts/geometry_continuation_study.py'
DEFAULT_CSV = ROOT / 'results/side_branches/ordering_principle/reference/geometry_continuation_study_first_test.csv'
args = sys.argv[1:]
if '--csv' not in args: args += ['--csv', str(DEFAULT_CSV)]
raise SystemExit(run_historical_script(SCRIPT, args))
