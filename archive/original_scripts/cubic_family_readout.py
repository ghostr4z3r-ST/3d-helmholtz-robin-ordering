import argparse
import importlib.util
import math
from pathlib import Path
from typing import Dict, Tuple, List

import numpy as np

# Load existing validated solver/diagnostic helpers
SRC = Path('/mnt/data/robin_3d_center_emergence.py')
spec = importlib.util.spec_from_file_location('rce', SRC)
rce = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rce)


def coeff_abs(sig: dict, key: str) -> float:
    return abs(sig['coeffs'].get(key, 0.0))


def average_abs_xyz(mesh, vec, locations: Dict[str, Tuple[float, float, float]], radius: float) -> float:
    vals = []
    for loc in locations.values():
        sig = rce.local_octant_signature(mesh, vec, loc, radius=radius)
        vals.append(coeff_abs(sig, 'xyz'))
    return float(np.mean(vals))


def readout(vals, vecs, mesh, radius: float):
    # center and six near-face centers, symmetric inside the cell
    eps = max(radius * 0.55, 0.12)
    eps = min(eps, 0.22)
    center = {'center': (0.5, 0.5, 0.5)}
    faces = {
        'face_x0': (eps, 0.5, 0.5),
        'face_x1': (1.0 - eps, 0.5, 0.5),
        'face_y0': (0.5, eps, 0.5),
        'face_y1': (0.5, 1.0 - eps, 0.5),
        'face_z0': (0.5, 0.5, eps),
        'face_z1': (0.5, 0.5, 1.0 - eps),
    }

    rows = []
    best_center = None
    best_faces = None
    for i in range(vecs.shape[1]):
        sig_c = rce.local_octant_signature(mesh, vecs[:, i], center['center'], radius=radius)
        xyz_center = coeff_abs(sig_c, 'xyz')
        xyz_faces = average_abs_xyz(mesh, vecs[:, i], faces, radius=radius)
        row = {
            'mode': i + 1,
            'lambda': float(vals[i]),
            'xyz_center': xyz_center,
            'xyz_faces_mean': xyz_faces,
            'center_top': sig_c['ranked'][:4],
        }
        rows.append(row)
        if best_center is None or xyz_center > best_center['xyz_center']:
            best_center = row
        if best_faces is None or xyz_faces > best_faces['xyz_faces_mean']:
            best_faces = row
    return rows, best_center, best_faces, center, faces


def print_summary(rows, best_center, best_faces):
    print('Mode summary (first modes):')
    for r in rows:
        top = ', '.join([f"{n}:{c:+.3f}" for n, c in r['center_top']])
        print(
            f"  mode {r['mode']:2d}: lambda={r['lambda']:.12f} | "
            f"|xyz|_center={r['xyz_center']:.6f} | "
            f"mean|xyz|_faces={r['xyz_faces_mean']:.6f} | center top -> {top}"
        )
    print('\nBest center-xyz mode:')
    print(
        f"  mode {best_center['mode']}, lambda={best_center['lambda']:.12f}, "
        f"|xyz|_center={best_center['xyz_center']:.6f}, "
        f"mean|xyz|_faces={best_center['xyz_faces_mean']:.6f}"
    )
    print('Best face-xyz mode:')
    print(
        f"  mode {best_faces['mode']}, lambda={best_faces['lambda']:.12f}, "
        f"|xyz|_center={best_faces['xyz_center']:.6f}, "
        f"mean|xyz|_faces={best_faces['xyz_faces_mean']:.6f}"
    )


def main():
    ap = argparse.ArgumentParser(description='Diagnostic readout for primitive-cubic vs bcc/fcc tendencies using existing cube solver.')
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--nx', type=int, default=10)
    ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--nz', type=int, default=10)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--radius', type=float, default=0.22)
    args = ap.parse_args()

    corners = rce.make_corners('cube')
    mesh = rce.build_hexahedral_mesh(args.nx, args.ny, args.nz, corners)
    A, M = rce.assemble(mesh, args.beta)
    vals, vecs = rce.solve_modes(A, M, args.modes)

    print(f'Cube diagnostic readout, beta={args.beta}')
    print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
    print(f'symmetry ||A-A^T||_F = {rce.symmetry_norm(A):.6e}')
    rows, best_center, best_faces, center, faces = readout(vals, vecs, mesh, args.radius)
    print_summary(rows, best_center, best_faces)


if __name__ == '__main__':
    main()
