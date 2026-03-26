from pathlib import Path
import sys
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/main_strand/nullmodel_blindtests/original_scripts/field_shuffle_supercell_null_tests_reconstructed.py'
DEFAULT_CSV = ROOT / 'results/main_strand/nullmodel_blindtests/reference/field_shuffle_supercell_null_tests.csv'
args = sys.argv[1:]
if '--csv' not in args: args += ['--csv', str(DEFAULT_CSV)]
raise SystemExit(run_historical_script(SCRIPT, args))
