
import argparse, importlib.util, math, csv
from pathlib import Path
import numpy as np

SRC = Path('/mnt/data/robin_3d_center_emergence.py')
spec = importlib.util.spec_from_file_location('rce', SRC)
rce = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rce)

def build_supercell_mesh(ncell: int, pts_per_cell: int):
    n = ncell * (pts_per_cell - 1) + 1
    corners = np.array([
        [0.0, 0.0, 0.0],
        [float(ncell), 0.0, 0.0],
        [float(ncell), float(ncell), 0.0],
        [0.0, float(ncell), 0.0],
        [0.0, 0.0, float(ncell)],
        [float(ncell), 0.0, float(ncell)],
        [float(ncell), float(ncell), float(ncell)],
        [0.0, float(ncell), float(ncell)],
    ], dtype=float)
    return rce.build_hexahedral_mesh(n, n, n, corners)

def subcell_centers(ncell: int):
    return [(i + 0.5, j + 0.5, k + 0.5) for k in range(ncell) for j in range(ncell) for i in range(ncell)]

def local_sig_physical(mesh, mode, center_phys, radius_phys):
    pts = mesh.points
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    L = maxs - mins
    logical_center = tuple(((np.array(center_phys) - mins) / L).tolist())
    logical_radius = radius_phys / float(L[0])
    return rce.local_octant_signature(mesh, mode, logical_center, radius=logical_radius)

def hadamard_basis_8():
    sigs = np.array([
        [-1, -1, -1], [ 1, -1, -1], [-1,  1, -1], [ 1,  1, -1],
        [-1, -1,  1], [ 1, -1,  1], [-1,  1,  1], [ 1,  1,  1],
    ], dtype=int)
    names = ['const', 'X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ']
    basis = []
    for name in names:
        if name == 'const':
            v = np.ones(8)
        elif name == 'X':
            v = sigs[:, 0]
        elif name == 'Y':
            v = sigs[:, 1]
        elif name == 'Z':
            v = sigs[:, 2]
        elif name == 'XY':
            v = sigs[:, 0] * sigs[:, 1]
        elif name == 'XZ':
            v = sigs[:, 0] * sigs[:, 2]
        elif name == 'YZ':
            v = sigs[:, 1] * sigs[:, 2]
        elif name == 'XYZ':
            v = sigs[:, 0] * sigs[:, 1] * sigs[:, 2]
        v = v.astype(float)
        v /= np.linalg.norm(v)
        basis.append(v)
    return names, np.vstack(basis)

def mode_metrics(mesh, vals, vecs, mode_index: int, ncell: int, radius_phys: float):
    centers = subcell_centers(ncell)
    coeffs_xyz = []
    axis_abs = []
    pair_abs = []
    for c in centers:
        sig = local_sig_physical(mesh, vecs[:, mode_index], c, radius_phys)
        coeffs = sig['coeffs']
        coeffs_xyz.append(coeffs['xyz'])
        axis_abs.append(max(abs(coeffs['x']), abs(coeffs['y']), abs(coeffs['z'])))
        pair_abs.append(max(abs(coeffs['xy']), abs(coeffs['xz']), abs(coeffs['yz'])))
    coeffs_xyz = np.asarray(coeffs_xyz, dtype=float)
    axis_abs = np.asarray(axis_abs, dtype=float)
    pair_abs = np.asarray(pair_abs, dtype=float)

    xyz_mean_abs = float(np.mean(np.abs(coeffs_xyz)))
    xyz_std_abs = float(np.std(np.abs(coeffs_xyz)))
    xyz_uniformity = float(1.0 - xyz_std_abs / (xyz_mean_abs + 1e-15))
    xyz_same_sign = float(abs(np.sum(coeffs_xyz)) / (np.sum(np.abs(coeffs_xyz)) + 1e-15))

    # default if no 8-cell hadamard readout is available
    xyz_const_pattern = xyz_same_sign
    best_nonconst_name = 'n/a'
    best_nonconst_abs = float('nan')

    # For 2x2x2 we can read the pattern class explicitly.
    if ncell == 2 and len(coeffs_xyz) == 8:
        names, basis = hadamard_basis_8()
        v = coeffs_xyz.copy()
        nrm = np.linalg.norm(v)
        vn = v / nrm if nrm > 0 else v
        patt_coeffs = basis @ vn
        patt = {name: float(patt_coeffs[i]) for i, name in enumerate(names)}
        xyz_const_pattern = abs(patt['const'])
        nonconst = {k: abs(vv) for k, vv in patt.items() if k != 'const'}
        best_nonconst_name = max(nonconst, key=nonconst.get)
        best_nonconst_abs = float(nonconst[best_nonconst_name])

    return {
        'mode': mode_index + 1,
        'lambda': float(vals[mode_index]),
        'xyz_mean_abs': xyz_mean_abs,
        'xyz_std_abs': xyz_std_abs,
        'xyz_uniformity': xyz_uniformity,
        'xyz_same_sign': xyz_same_sign,
        'xyz_const_pattern': xyz_const_pattern,
        'xyz_best_nonconst_pattern': best_nonconst_name,
        'xyz_best_nonconst_abs': best_nonconst_abs,
        'axis_mean_abs': float(np.mean(axis_abs)),
        'pair_mean_abs': float(np.mean(pair_abs)),
        'dominance_over_axis': float(xyz_mean_abs / (float(np.mean(axis_abs)) + 1e-15)),
        'dominance_over_pair': float(xyz_mean_abs / (float(np.mean(pair_abs)) + 1e-15)),
    }

def score_row(r, pattern_weight='const'):
    if pattern_weight == 'const':
        pattern = r['xyz_const_pattern']
    else:
        alt = r['xyz_best_nonconst_abs']
        if isinstance(alt, float) and np.isnan(alt):
            pattern = r['xyz_const_pattern']
        else:
            pattern = max(r['xyz_const_pattern'], alt)
    return r['xyz_mean_abs'] * pattern * max(r['xyz_uniformity'], 0.0)

def analyze(beta: float, ncell: int, pts_per_cell: int, modes: int, radius: float, skip_ground=True, pattern_weight='const'):
    mesh = build_supercell_mesh(ncell, pts_per_cell)
    A, M = rce.assemble(mesh, beta)
    vals, vecs = rce.solve_modes(A, M, modes)
    start = 1 if skip_ground else 0
    rows = [mode_metrics(mesh, vals, vecs, i, ncell, radius) for i in range(start, vecs.shape[1])]
    for r in rows:
        r['score'] = score_row(r, pattern_weight=pattern_weight)
        r['beta'] = beta
        r['ncell'] = ncell
        r['pts_per_cell'] = pts_per_cell
    best = max(rows, key=lambda r: r['score'])
    return mesh, vals, vecs, rows, best, rce.symmetry_norm(A)

def print_summary(beta, rows, best, max_rows=12):
    print(f'=== primitive cubic supercell: beta={beta} ===')
    print('mode | lambda | mean|xyz| | uniform | same-sign | const-pattern | best-cell-pattern | axis | pair | score')
    for r in rows[:max_rows]:
        alt = r['xyz_best_nonconst_abs']
        alt_text = f"{r['xyz_best_nonconst_pattern']}:{alt:.3f}" if not (isinstance(alt,float) and np.isnan(alt)) else 'n/a'
        print(f"{r['mode']:4d} | {r['lambda']:.9f} | {r['xyz_mean_abs']:.6f} | {r['xyz_uniformity']:.3f} | "
              f"{r['xyz_same_sign']:.3f} | {r['xyz_const_pattern']:.3f} | {alt_text} | "
              f"{r['axis_mean_abs']:.3f} | {r['pair_mean_abs']:.3f} | {r['score']:.6f}")
    print('best nontrivial ordered xyz mode:')
    print(best)

def main():
    ap = argparse.ArgumentParser(description='Hardened primitive cubic supercell xyz-order diagnostics.')
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5,10')
    ap.add_argument('--ncell', type=int, default=2)
    ap.add_argument('--pts-per-cell', type=int, default=8)
    ap.add_argument('--modes', type=int, default=24)
    ap.add_argument('--radius', type=float, default=0.35)
    ap.add_argument('--include-ground', action='store_true')
    ap.add_argument('--pattern-weight', choices=['const', 'any'], default='const',
                    help='const: prefer same-sign translational order; any: allow alternating cell patterns too.')
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = [float(x) for x in args.betas.split(',') if x.strip()]
    allbest = []
    for beta in betas:
        mesh, vals, vecs, rows, best, sym = analyze(beta, args.ncell, args.pts_per_cell, args.modes, args.radius,
                                                    skip_ground=not args.include_ground,
                                                    pattern_weight=args.pattern_weight)
        print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
        print(f'symmetry ||A-A^T||_F = {sym:.6e}')
        print_summary(beta, rows, best, max_rows=min(12, len(rows)))
        print()
        allbest.append(best)

    if args.csv:
        fields = ['beta','mode','lambda','xyz_mean_abs','xyz_std_abs','xyz_uniformity','xyz_same_sign',
                  'xyz_const_pattern','xyz_best_nonconst_pattern','xyz_best_nonconst_abs',
                  'axis_mean_abs','pair_mean_abs','dominance_over_axis','dominance_over_pair','score',
                  'ncell','pts_per_cell']
        with open(args.csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in allbest:
                w.writerow({k: r[k] for k in fields})

if __name__ == '__main__':
    main()
