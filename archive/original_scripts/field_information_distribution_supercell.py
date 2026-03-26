from __future__ import annotations

import csv
import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np

# Reuse the existing supercell machinery and the already established mode selections
SRC_RCE = Path('/mnt/data/robin_3d_center_emergence.py')
SRC_PCS = Path('/mnt/data/primitive_cubic_supercell_ordered.py')
PHYS_CSV = Path('/mnt/data/cubic_family_physical_readout_results.csv')

spec_rce = importlib.util.spec_from_file_location('rce', SRC_RCE)
rce = importlib.util.module_from_spec(spec_rce)
spec_rce.loader.exec_module(rce)

spec_pcs = importlib.util.spec_from_file_location('pcs', SRC_PCS)
pcs = importlib.util.module_from_spec(spec_pcs)
spec_pcs.loader.exec_module(pcs)


@dataclass
class Row:
    beta: float
    family: str
    mode_index: int
    eigenvalue: float
    f_eff: float
    entropy_norm: float
    anisotropy_ratio: float
    eig1: float
    eig2: float
    eig3: float
    boundary_ratio: float
    note: str


def row_sum_weights(M) -> np.ndarray:
    return np.asarray(M @ np.ones(M.shape[0])).ravel()


def load_mode_selections() -> Dict[float, Dict[str, int]]:
    selections: Dict[float, Dict[str, int]] = {}
    with PHYS_CSV.open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            beta = float(r['beta'])
            if beta in (0.0, 1.0, 5.0):
                selections[beta] = {
                    'primitive': int(r['primitive_mode']) - 1,
                    'fcc': int(r['fcc_mode']) - 1,
                }
    return selections


def field_distribution_metrics(mesh, modevec: np.ndarray, weights: np.ndarray):
    # normalize in the lumped-mass measure
    rho = np.asarray(modevec, dtype=float) ** 2
    rho = rho / (np.sum(weights * rho) + 1e-15)

    total_volume = float(np.sum(weights))
    ipr = float(np.sum(weights * rho * rho))
    f_eff = float(1.0 / (ipr * total_volume))

    p = weights * rho
    p = p / p.sum()
    eps = 1e-300
    entropy_norm = float(-np.sum(p * np.log(p + eps)) / np.log(len(p)))

    pts = mesh.points
    center = np.array([
        np.sum(weights * rho * pts[:, 0]),
        np.sum(weights * rho * pts[:, 1]),
        np.sum(weights * rho * pts[:, 2]),
    ])
    d = pts - center[None, :]
    M2 = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            M2[a, b] = np.sum(weights * rho * d[:, a] * d[:, b])
    eigs = np.sort(np.linalg.eigvalsh(M2))[::-1]
    anis = float(eigs[-1] / eigs[0]) if eigs[0] > 0 else 0.0

    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    on_boundary = np.any(np.isclose(pts, mins) | np.isclose(pts, maxs), axis=1)
    boundary_ratio = float(np.sum(weights[on_boundary] * rho[on_boundary]))
    return f_eff, entropy_norm, anis, eigs, boundary_ratio


def run_test(ncell: int = 3, pts_per_cell: int = 6, betas: List[float] | None = None, modes: int = 20) -> List[Row]:
    if betas is None:
        betas = [0.0, 1.0, 5.0]
    selections = load_mode_selections()
    mesh = pcs.build_supercell_mesh(ncell, pts_per_cell)

    rows: List[Row] = []
    for beta in betas:
        A, M = rce.assemble(mesh, beta)
        vals, vecs = rce.solve_modes(A, M, modes)
        weights = row_sum_weights(M)
        prim_mode = selections[beta]['primitive']
        fcc_mode = selections[beta]['fcc']
        for family, mi in [('primitive_q', prim_mode), ('fcc_face', fcc_mode)]:
            f_eff, H, anis, eigs, B = field_distribution_metrics(mesh, vecs[:, mi], weights)
            note = 'same_underlying_mode_as_other_family' if prim_mode == fcc_mode else 'distinct_mode'
            rows.append(Row(beta, family, mi + 1, float(vals[mi]), f_eff, H, anis, float(eigs[0]), float(eigs[1]), float(eigs[2]), B, note))
    return rows


def write_csv(rows: List[Row], path: Path) -> None:
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'beta', 'family', 'mode_index', 'eigenvalue', 'f_eff', 'entropy_norm',
            'anisotropy_ratio', 'eig1', 'eig2', 'eig3', 'boundary_ratio', 'note'
        ])
        for r in rows:
            writer.writerow([
                r.beta, r.family, r.mode_index, r.eigenvalue, r.f_eff, r.entropy_norm,
                r.anisotropy_ratio, r.eig1, r.eig2, r.eig3, r.boundary_ratio, r.note
            ])


def main() -> None:
    out = Path('/mnt/data/field_information_distribution_supercell_test.csv')
    rows = run_test()
    write_csv(rows, out)
    print(f'wrote {out}')
    for r in rows:
        print(r)


if __name__ == '__main__':
    main()
