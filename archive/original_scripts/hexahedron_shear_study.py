import argparse
import csv
from pathlib import Path
from typing import List

from robin_3d_center_emergence import (
    make_corners,
    build_hexahedral_mesh,
    assemble,
    solve_modes,
    symmetry_norm,
    first_xyz_mode,
    local_octant_signature,
)


def run_case(nx:int, ny:int, nz:int, beta:float, shear_xy:float, modes:int, radius:float,
             Lx:float=1.0, Ly:float=1.0, Lz:float=1.0):
    corners = make_corners('hexahedron', Lx, Ly, Lz, shear_xy=shear_xy, shear_xz=0.0, shear_yz=0.0)
    mesh = build_hexahedral_mesh(nx, ny, nz, corners)
    A, M = assemble(mesh, beta)
    vals, vecs = solve_modes(A, M, modes)

    locs = {
        'center': (0.5, 0.5, 0.5),
        'face_x0': (0.12, 0.5, 0.5),
        'edge_x0y0': (0.12, 0.12, 0.5),
        'corner_near': (0.12, 0.12, 0.12),
    }
    best_idx, xyz_amp, center_sig = first_xyz_mode(mesh, vecs, logical_center=locs['center'], radius=radius)

    face_sig = local_octant_signature(mesh, vecs[:, best_idx], locs['face_x0'], radius=radius)
    edge_sig = local_octant_signature(mesh, vecs[:, best_idx], locs['edge_x0y0'], radius=radius)
    corner_sig = local_octant_signature(mesh, vecs[:, best_idx], locs['corner_near'], radius=radius)

    return {
        'shear_xy': shear_xy,
        'beta': beta,
        'symmetry_norm': symmetry_norm(A),
        'best_mode': best_idx + 1,
        'lambda_best': float(vals[best_idx]),
        'xyz_center': abs(center_sig['coeffs']['xyz']),
        'xyz_face': abs(face_sig['coeffs']['xyz']),
        'xyz_edge': abs(edge_sig['coeffs']['xyz']),
        'xyz_corner': abs(corner_sig['coeffs']['xyz']),
        'center_top1': center_sig['ranked'][0][0],
        'center_top1_coeff': center_sig['ranked'][0][1],
        'center_top2': center_sig['ranked'][1][0],
        'center_top2_coeff': center_sig['ranked'][1][1],
    }


def parse_floats(text: str) -> List[float]:
    return [float(x) for x in text.split(',') if x.strip()]


def main():
    ap = argparse.ArgumentParser(description='Shear-family study for hexahedral Robin/Helmholtz center-emergence diagnostics.')
    ap.add_argument('--nx', type=int, default=10)
    ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--nz', type=int, default=10)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--radius', type=float, default=0.22)
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5')
    ap.add_argument('--shears', type=str, default='0.0,0.05,0.10,0.15,0.20,0.30')
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = parse_floats(args.betas)
    shears = parse_floats(args.shears)

    rows = []
    for beta in betas:
        print(f'### beta={beta}')
        for shear in shears:
            row = run_case(args.nx, args.ny, args.nz, beta, shear, args.modes, args.radius)
            rows.append(row)
            print(
                f"shear={shear:0.3f} | best_mode={row['best_mode']:2d} | "
                f"lambda={row['lambda_best']:.6f} | xyz_center={row['xyz_center']:.3f} | "
                f"xyz_face={row['xyz_face']:.3f} | xyz_edge={row['xyz_edge']:.3f} | xyz_corner={row['xyz_corner']:.3f} | "
                f"top1={row['center_top1']}({row['center_top1_coeff']:+.3f}) top2={row['center_top2']}({row['center_top2_coeff']:+.3f})"
            )
        print()

    if args.csv:
        path = Path(args.csv)
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f'Wrote CSV: {path}')


if __name__ == '__main__':
    main()
