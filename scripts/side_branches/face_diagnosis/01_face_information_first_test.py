from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT / 'src') not in sys.path:
    sys.path.insert(0, str(ROOT / 'src'))

from paper1.face_information import main

if __name__ == '__main__':
    # Default output for the side-branch subfolder when no explicit --csv is given.
    if '--csv' not in sys.argv:
        sys.argv.extend([
            '--csv',
            str(ROOT / 'results/side_branches/face_diagnosis/reference/face_information_first_test_rerun.csv')
        ])
    main()
