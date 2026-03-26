from __future__ import annotations

"""Best-effort reconstruction of the missing historical script
`field_shuffle_supercell_null_tests.py`.

This reconstruction is based on the surviving neighboring scripts and result CSVs:
- field_shuffle_null_tests.py (cube field-shuffle pattern)
- face_information_sidebranch_extended.py
- face_patch_scan_supercell.py
- primitive_cubic_supercell_ordered.py
- fcc_supercell_readout.py
- field_shuffle_supercell_null_tests.csv
- supercell_qdominated_field_shuffle.csv

It aims to recover the historical *logic* of the test:
- select primitive finite-q and fcc-like top modes on a 3x3x3 supercell
- shuffle field values nodewise while preserving the amplitude distribution
- recompute outer-face and volume readouts on the shuffled fields
- compare real vs. field-shuffled classification accuracies
- compare key q/order metrics against the shuffled null distribution

Because the original source file is missing and the surviving stub is empty, this
script should be treated as a transparent reconstruction rather than a guaranteed
byte-identical restoration.
"""

import argparse
import csv
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import numpy as np


ROOT = Path('/mnt/data')
FACES = ['x0', 'x1', 'y0', 'y1', 'z0', 'z1']


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


fise = _load('fise_recon', ROOT / 'face_information_sidebranch_extended.py')
fpss = _load('fpss_recon', ROOT / 'face_patch_scan_supercell.py')
pcs = _load('pcs_recon', ROOT / 'primitive_cubic_supercell_ordered.py')
fccr = _load('fccr_recon', ROOT / 'fcc_supercell_readout.py')


@dataclass
class SelectedMode:
    family: str
    beta: float
    mesh: object
    mode_index: int
    modevec: np.ndarray
    row: Dict[str, float]


def rng_from_seed(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)


def shuffle_mode(modevec: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    return np.asarray(modevec, float)[rng.permutation(len(modevec))]


def primitive_row_from_mode(mesh, modevec: np.ndarray, ncell: int = 3, radius: float = 0.28) -> Dict[str, float]:
    vals = np.array([np.nan], dtype=float)
    vecs = np.asarray(modevec, float).reshape(-1, 1)
    row = pcs.mode_metrics(mesh, vals, vecs, 0, ncell, radius)
    row['score'] = pcs.score_row(row)
    return row


def fcc_row_from_mode(mesh, modevec: np.ndarray, ncell: int = 3, radius: float = 0.28) -> Dict[str, float]:
    vals = np.array([np.nan], dtype=float)
    vecs = np.asarray(modevec, float).reshape(-1, 1)
    return fccr.mode_metrics(mesh, vals, vecs, 0, ncell, radius)


def select_modes(
    betas: Sequence[float] = (0.0, 1.0, 5.0),
    ncell: int = 3,
    pts_per_cell: int = 6,
    modes: int = 20,
    radius: float = 0.28,
    topk: int = 3,
) -> List[SelectedMode]:
    selected: List[SelectedMode] = []
    for beta in betas:
        meshp, valsp, vecsp, rowsp, bestp, symp = pcs.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_prim = sorted(rowsp, key=lambda r: r['score'], reverse=True)[:topk]
        for row in top_prim:
            mi = row['mode'] - 1
            selected.append(SelectedMode('primitive_q', beta, meshp, mi, vecsp[:, mi].copy(), dict(row)))

        meshf, valsf, vecsf, rowsf, bestc, bestf, symf = fccr.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_fcc = sorted(rowsf, key=lambda r: r['fcc_score'], reverse=True)[:topk]
        for row in top_fcc:
            mi = row['mode'] - 1
            selected.append(SelectedMode('fcc_face', beta, meshf, mi, vecsf[:, mi].copy(), dict(row)))
    return selected


def build_classification_samples(
    selected: Sequence[SelectedMode],
    shuffle: bool,
    rng: np.random.Generator,
    face_radius: float = 0.30,
    face_thickness: float = 0.08,
    radius_phys: float = 0.28,
) -> Dict[str, List[Dict[str, np.ndarray]]]:
    """Build the three supercell classification tasks used in the archived CSV.

    face_radius=0.30 is chosen because it matches the archived supercell face
    accuracies more closely than the cube defaults while keeping the surviving
    volume and q-metric logic unchanged.
    """
    raw_1f: List[Dict[str, np.ndarray]] = []
    raw_6f: List[Dict[str, np.ndarray]] = []
    vols: List[Dict[str, np.ndarray]] = []

    for item in selected:
        mode = shuffle_mode(item.modevec, rng) if shuffle else item.modevec

        fsigs = fpss.supercell_outer_face_sigs(item.mesh, mode, radius=face_radius, thickness=face_thickness)
        for f in FACES:
            raw_1f.append({'label': item.family, 'feature': fise.improved_face_vector(fsigs[f])})
        raw_6f.append({'label': item.family, 'feature': fise.improved_sixface_vector(fsigs)})

        if item.family == 'primitive_q':
            row = primitive_row_from_mode(item.mesh, mode, radius=radius_phys)
            vols.append({'label': item.family, 'feature': fpss.primitive_volume_feature(row)})
        else:
            row = fcc_row_from_mode(item.mesh, mode, radius=radius_phys)
            vols.append({'label': item.family, 'feature': fpss.fcc_volume_feature(row)})

    return {
        'supercell_raw_1F': raw_1f,
        'supercell_raw_6F': raw_6f,
        'supercell_volume': vols,
    }


def classification_accuracy(samples: List[Dict[str, np.ndarray]]) -> float:
    return float(fise.nearest_centroid_cv(samples))


def q_metric_means(
    selected: Sequence[SelectedMode],
    shuffle: bool,
    rng: np.random.Generator,
    radius_phys: float = 0.28,
) -> Dict[str, float]:
    prim_rows: List[Dict[str, float]] = []
    fcc_rows: List[Dict[str, float]] = []
    for item in selected:
        mode = shuffle_mode(item.modevec, rng) if shuffle else item.modevec
        if item.family == 'primitive_q':
            prim_rows.append(primitive_row_from_mode(item.mesh, mode, radius=radius_phys))
        else:
            fcc_rows.append(fcc_row_from_mode(item.mesh, mode, radius=radius_phys))

    out = {
        'primitive_q:xyz_mean_abs': float(np.mean([r['xyz_mean_abs'] for r in prim_rows])),
        'primitive_q:xyz_uniformity': float(np.mean([r['xyz_uniformity'] for r in prim_rows])),
        'primitive_q:best_q_amp': float(np.mean([r['best_q_amp'] for r in prim_rows])),
        'primitive_q:score': float(np.mean([r['score'] for r in prim_rows])),
        'fcc:fcc_face_mean': float(np.mean([r['fcc_face_mean'] for r in fcc_rows])),
        'fcc:fcc_face_uniformity': float(np.mean([r['fcc_face_uniformity'] for r in fcc_rows])),
        'fcc:fcc_face_qamp': float(np.mean([r['fcc_face_qamp'] for r in fcc_rows])),
        'fcc:fcc_score': float(np.mean([r['fcc_score'] for r in fcc_rows])),
        'fcc:center_score': float(np.mean([r['center_score'] for r in fcc_rows])),
    }
    return out


def summarize_null(real: float, shuffled: Sequence[float]) -> Dict[str, float]:
    arr = np.asarray(shuffled, dtype=float)
    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0
    z = float((real - mean) / (std + 1e-15))
    p_ge = float((np.sum(arr >= real) + 1) / (len(arr) + 1))
    return {
        'real': float(real),
        'shuffle_mean': mean,
        'shuffle_std': std,
        'z_score': z,
        'p_ge': p_ge,
        'n_perm': int(len(arr)),
    }


def run_test(
    n_perm: int = 20,
    seed: int = 12345,
    betas: Sequence[float] = (0.0, 1.0, 5.0),
    ncell: int = 3,
    pts_per_cell: int = 6,
    modes: int = 20,
    radius: float = 0.28,
    topk: int = 3,
    face_radius: float = 0.30,
    face_thickness: float = 0.08,
) -> List[Dict[str, float]]:
    selected = select_modes(betas=betas, ncell=ncell, pts_per_cell=pts_per_cell, modes=modes, radius=radius, topk=topk)

    real_rng = rng_from_seed(seed)
    real_class = build_classification_samples(selected, shuffle=False, rng=real_rng, face_radius=face_radius, face_thickness=face_thickness, radius_phys=radius)
    real_class_acc = {k: classification_accuracy(v) for k, v in real_class.items()}
    real_q = q_metric_means(selected, shuffle=False, rng=real_rng, radius_phys=radius)

    shuf_class: Dict[str, List[float]] = {k: [] for k in real_class_acc}
    shuf_q: Dict[str, List[float]] = {k: [] for k in real_q}

    # Independent permutations per replicate, but shared across the metrics within one replicate.
    for rep in range(n_perm):
        rng = rng_from_seed(seed + rep + 1)
        cls = build_classification_samples(selected, shuffle=True, rng=rng, face_radius=face_radius, face_thickness=face_thickness, radius_phys=radius)
        acc = {k: classification_accuracy(v) for k, v in cls.items()}
        for k, v in acc.items():
            shuf_class[k].append(v)

        rng_q = rng_from_seed(seed + 10000 + rep + 1)
        qvals = q_metric_means(selected, shuffle=True, rng=rng_q, radius_phys=radius)
        for k, v in qvals.items():
            shuf_q[k].append(v)

    rows: List[Dict[str, float]] = []
    for task in ['supercell_raw_1F', 'supercell_raw_6F', 'supercell_volume']:
        row = {'section': 'classification', 'task': task}
        row.update(summarize_null(real_class_acc[task], shuf_class[task]))
        rows.append(row)

    for task in [
        'primitive_q:xyz_mean_abs',
        'primitive_q:xyz_uniformity',
        'primitive_q:best_q_amp',
        'primitive_q:score',
        'fcc:fcc_face_mean',
        'fcc:fcc_face_uniformity',
        'fcc:fcc_face_qamp',
        'fcc:fcc_score',
        'fcc:center_score',
    ]:
        row = {'section': 'q_metrics', 'task': task}
        row.update(summarize_null(real_q[task], shuf_q[task]))
        rows.append(row)

    return rows


def write_csv(rows: Sequence[Dict[str, float]], path: Path) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        fieldnames = ['section', 'task', 'real', 'shuffle_mean', 'shuffle_std', 'z_score', 'p_ge', 'n_perm']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser(description='Best-effort reconstruction of the missing supercell field-shuffle null test.')
    ap.add_argument('--n-perm', type=int, default=20)
    ap.add_argument('--seed', type=int, default=12345)
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--radius', type=float, default=0.28)
    ap.add_argument('--topk', type=int, default=3)
    ap.add_argument('--face-radius', type=float, default=0.30)
    ap.add_argument('--face-thickness', type=float, default=0.08)
    ap.add_argument('--csv', type=str, default='/mnt/data/field_shuffle_supercell_null_tests_reconstructed.csv')
    args = ap.parse_args()

    rows = run_test(
        n_perm=args.n_perm,
        seed=args.seed,
        ncell=args.ncell,
        pts_per_cell=args.pts_per_cell,
        modes=args.modes,
        radius=args.radius,
        topk=args.topk,
        face_radius=args.face_radius,
        face_thickness=args.face_thickness,
    )
    out = Path(args.csv)
    write_csv(rows, out)
    print(out)
    for r in rows:
        print(r)


if __name__ == '__main__':
    main()
