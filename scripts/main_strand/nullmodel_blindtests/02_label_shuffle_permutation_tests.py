from pathlib import Path
import sys
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/main_strand/nullmodel_blindtests/original_scripts/label_shuffle_permutation_tests.py'
DEFAULT_CSV = ROOT / 'results/main_strand/nullmodel_blindtests/reference/label_shuffle_permutation_tests.csv'
args = sys.argv[1:]
if '--csv' not in args: args += ['--csv', str(DEFAULT_CSV)]
raise SystemExit(run_historical_script(SCRIPT, args))
