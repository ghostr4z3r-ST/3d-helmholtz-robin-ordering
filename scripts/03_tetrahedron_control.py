from __future__ import annotations

import argparse
import csv
from pathlib import Path

from paper1.tetra_fem import regular_tetra_vertices, solve, local_octant_signature


def main() -> None:
    ap = argparse.ArgumentParser(description='Tetrahedron control case for center-signature comparison.')
    ap.add_argument('--n', type=int, default=10)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    verts = regular_tetra_vertices()
    center = verts.mean(axis=0)
    res = solve(verts, n=args.n, beta=args.beta, modes=args.modes)
    rows = []
    best_mode = None
    best_xyz = -1.0
    for i, lam in enumerate(res.eigenvalues):
        coeffs = local_octant_signature(res.points, res.eigenvectors[:, i], center)
        row = {'mode': i + 1, 'lambda': float(lam), **{k: float(v) for k, v in coeffs.items()}}
        rows.append(row)
        if abs(coeffs['xyz']) > best_xyz:
            best_xyz = abs(coeffs['xyz'])
            best_mode = i + 1
    print(f'best_xyz_mode={best_mode}, abs_c_xyz={best_xyz:.6f}')
    for row in rows[: min(8, len(rows))]:
        print(row)

    if args.csv:
        path = Path(args.csv)
        with path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)


if __name__ == '__main__':
    main()
