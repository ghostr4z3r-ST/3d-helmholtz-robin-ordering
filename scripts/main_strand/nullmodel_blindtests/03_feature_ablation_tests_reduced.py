from pathlib import Path
from paper1.historical_runner import run_historical_script
ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / 'archive/main_strand/nullmodel_blindtests/original_scripts/feature_ablation_tests_reduced.py'
CAPTURE = {'feature_ablation_tests_reduced.csv': ROOT / 'results/main_strand/nullmodel_blindtests/reference/feature_ablation_tests_reduced.csv', 'feature_ablation_summary_reduced.csv': ROOT / 'results/main_strand/nullmodel_blindtests/reference/feature_ablation_summary_reduced.csv'}
raise SystemExit(run_historical_script(SCRIPT, [], capture_outputs=CAPTURE))
