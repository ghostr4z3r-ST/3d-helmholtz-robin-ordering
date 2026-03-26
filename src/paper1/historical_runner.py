from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Sequence


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _stage_known_files(stage_dir: Path) -> None:
    root = repo_root()
    exts = {'.py', '.csv', '.png', '.txt', '.md'}
    seen: set[str] = set()
    for path in root.rglob('*'):
        if not path.is_file() or path.suffix.lower() not in exts:
            continue
        name = path.name
        if name in seen:
            continue
        dst = stage_dir / name
        try:
            os.symlink(path, dst)
        except OSError:
            shutil.copy2(path, dst)
        seen.add(name)


def run_historical_script(script_path: Path, argv: Sequence[str] | None = None,
                          capture_outputs: Mapping[str, Path] | None = None) -> int:
    argv = list(argv or [])
    capture_outputs = dict(capture_outputs or {})
    script_path = Path(script_path).resolve()
    if not script_path.exists():
        raise FileNotFoundError(script_path)

    with tempfile.TemporaryDirectory(prefix='paper1_hist_') as td:
        stage = Path(td)
        _stage_known_files(stage)
        text = script_path.read_text(encoding='utf-8', errors='ignore')
        text = text.replace('/mnt/data/', str(stage) + '/')
        staged_script = stage / script_path.name
        staged_script.write_text(text, encoding='utf-8')
        completed = subprocess.run([sys.executable, str(staged_script), *argv], cwd=stage)
        if completed.returncode != 0:
            return completed.returncode
        for basename, target in capture_outputs.items():
            src = stage / basename
            if src.exists():
                target = Path(target)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, target)
        return completed.returncode
