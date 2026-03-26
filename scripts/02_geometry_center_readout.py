from __future__ import annotations

import argparse
import csv
from pathlib import Path

from paper1.hexa_fem import make_corners, build_hexahedral_mesh, assemble, solve_modes, symmetry_norm
from paper1.octant import first_xyz_mode, local_octant_signature


def coeff_abs(sig: dict, key: str) -> float:
    return abs(sig['coeffs'].get(key, 0.0))


def run_case(shape: str, beta: float, nx: int, ny: int, nz: int, modes: int, radius: float,
             Lx: float = 1.0, Ly: float = 1.0, Lz: float = 1.0,
             shear_xy: float = 0.0, shear_xz: float = 0.0, shear_yz: float = 0.0) -> dict:
    corners = make_corners(shape, Lx=Lx, Ly=Ly, Lz=Lz, shear_xy=shear_xy, shear_xz=shear_xz, shear_yz=shear_yz)
    mesh = build_hexahedral_mesh(nx, ny, nz, corners)
    A, M, _ = assemble(mesh, beta)
    vals, vecs = solve_modes(A, M, modes)
    best_idx, best_val, best_sig = first_xyz_mode(mesh, vecs, logical_center=(0.5, 0.5, 0.5), radius=radius)
    face_sig = local_octant_signature(mesh, vecs[:, best_idx], logical_center=(0.12, 0.5, 0.5), radius=radius)
    return {
        'shape': shape,
        'beta': beta,
        'best_mode': best_idx + 1,
        'lambda_best': float(vals[best_idx]),
        'xyz_center': float(best_val),
        'xyz_face': coeff_abs(face_sig, 'xyz'),
        'center_top1': best_sig['ranked'][0][0],
        'center_top1_coeff': float(best_sig['ranked'][0][1]),
        'center_top2': best_sig['ranked'][1][0],
        'center_top2_coeff': float(best_sig['ranked'][1][1]),
        'symmetry_norm': symmetry_norm(A),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description='Compare center-emergence diagnostics across cube/quader/hexahedron.')
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--nx', type=int, default=10)
    ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--nz', type=int, default=10)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--radius', type=float, default=0.22)
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    cases = [
        ('cube', dict()),
        ('quader', dict(Lx=1.0, Ly=1.2, Lz=0.8)),
        ('hexahedron', dict(Lx=1.0, Ly=1.0, Lz=1.0, shear_xy=0.15)),
    ]
    rows = []
    for shape, kw in cases:
        row = run_case(shape, args.beta, args.nx, args.ny, args.nz, args.modes, args.radius, **kw)
        rows.append(row)
        print(row)

    if args.csv:
        path = Path(args.csv)
        with path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)


if __name__ == '__main__':
    main()
