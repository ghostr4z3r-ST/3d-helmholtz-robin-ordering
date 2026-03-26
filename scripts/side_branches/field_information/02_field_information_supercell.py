from pathlib import Path
import sys
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/side_branches/field_information/original_scripts/field_information_distribution_supercell.py'
DEFAULT_CSV = ROOT / 'results/side_branches/field_information/reference/field_information_distribution_supercell_test.csv'
args = sys.argv[1:]
if '--csv' not in args: args += ['--csv', str(DEFAULT_CSV)]
raise SystemExit(run_historical_script(SCRIPT, args))
