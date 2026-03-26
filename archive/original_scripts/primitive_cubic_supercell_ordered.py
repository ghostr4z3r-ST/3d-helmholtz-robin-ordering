import argparse, importlib.util, csv, math
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


def lattice_fourier_diagnostics(coeffs_xyz: np.ndarray, ncell: int):
    # coeffs_xyz shaped (ncell,ncell,ncell) indexed [k,j,i]
    total_abs = float(np.sum(np.abs(coeffs_xyz))) + 1e-15
    best_amp = -1.0
    best_q = (0, 0, 0)
    q_amps = {}
    for mx in range(ncell):
        for my in range(ncell):
            for mz in range(ncell):
                s = 0.0 + 0.0j
                for k in range(ncell):
                    for j in range(ncell):
                        for i in range(ncell):
                            phase = -2j * math.pi * (mx * i + my * j + mz * k) / ncell
                            s += coeffs_xyz[k, j, i] * np.exp(phase)
                amp = abs(s) / total_abs
                q = (mx, my, mz)
                q_amps[q] = float(amp)
                if amp > best_amp:
                    best_amp = float(amp)
                    best_q = q
    return best_q, best_amp, q_amps


def neighbor_corr(coeffs_xyz: np.ndarray, axis: int):
    # axis: 0->x, 1->y, 2->z over periodic neighbor pairs
    total = 0.0
    denom = 0.0
    nz = coeffs_xyz.shape[0]
    ny = coeffs_xyz.shape[1]
    nx = coeffs_xyz.shape[2]
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                a = coeffs_xyz[k, j, i]
                if axis == 0:
                    b = coeffs_xyz[k, j, (i + 1) % nx]
                elif axis == 1:
                    b = coeffs_xyz[k, (j + 1) % ny, i]
                else:
                    b = coeffs_xyz[(k + 1) % nz, j, i]
                total += a * b
                denom += abs(a * b)
    if denom < 1e-15:
        return 0.0
    return float(total / denom)


def classify_q(q, ncell):
    if q == (0, 0, 0):
        return 'const'
    # for ncell=2 this recovers X/Y/Z... logic; for ncell=3 use q tuple labels.
    names = []
    labels = ['X', 'Y', 'Z']
    for comp, lab in zip(q, labels):
        if comp % ncell != 0:
            names.append(f'{lab}{comp}')
    return '+'.join(names) if names else 'const'


def mode_metrics(mesh, vals, vecs, mode_index: int, ncell: int, radius_phys: float):
    centers = subcell_centers(ncell)
    coeffs_xyz = []
    coeffs_axis = []
    coeffs_pair = []
    for c in centers:
        sig = local_sig_physical(mesh, vecs[:, mode_index], c, radius_phys)
        coeffs = sig['coeffs']
        coeffs_xyz.append(coeffs['xyz'])
        coeffs_axis.append(max(abs(coeffs['x']), abs(coeffs['y']), abs(coeffs['z'])))
        coeffs_pair.append(max(abs(coeffs['xy']), abs(coeffs['xz']), abs(coeffs['yz'])))
    coeffs_xyz = np.asarray(coeffs_xyz, dtype=float)
    coeffs_axis = np.asarray(coeffs_axis, dtype=float)
    coeffs_pair = np.asarray(coeffs_pair, dtype=float)
    cube = coeffs_xyz.reshape((ncell, ncell, ncell))

    xyz_mean_abs = float(np.mean(np.abs(coeffs_xyz)))
    xyz_std_abs = float(np.std(np.abs(coeffs_xyz)))
    xyz_uniformity = float(max(0.0, 1.0 - xyz_std_abs / (xyz_mean_abs + 1e-15)))
    xyz_same_sign = float(abs(np.sum(coeffs_xyz)) / (np.sum(np.abs(coeffs_xyz)) + 1e-15))

    best_q, best_q_amp, _ = lattice_fourier_diagnostics(cube, ncell)
    corr_x = neighbor_corr(cube, 0)
    corr_y = neighbor_corr(cube, 1)
    corr_z = neighbor_corr(cube, 2)
    corr_mean = float((corr_x + corr_y + corr_z) / 3.0)

    return {
        'mode': mode_index + 1,
        'lambda': float(vals[mode_index]),
        'xyz_mean_abs': xyz_mean_abs,
        'xyz_std_abs': xyz_std_abs,
        'xyz_uniformity': xyz_uniformity,
        'xyz_same_sign': xyz_same_sign,
        'best_q': best_q,
        'best_q_label': classify_q(best_q, ncell),
        'best_q_amp': best_q_amp,
        'corr_x': corr_x,
        'corr_y': corr_y,
        'corr_z': corr_z,
        'corr_mean': corr_mean,
        'axis_mean_abs': float(np.mean(coeffs_axis)),
        'pair_mean_abs': float(np.mean(coeffs_pair)),
        'dominance_over_axis': float(xyz_mean_abs / (float(np.mean(coeffs_axis)) + 1e-15)),
        'dominance_over_pair': float(xyz_mean_abs / (float(np.mean(coeffs_pair)) + 1e-15)),
    }


def score_row(r):
    # seek strong, uniform, globally ordered xyz; reward dominant lattice-order pattern
    return (
        r['xyz_mean_abs']
        * r['xyz_uniformity']
        * max(r['best_q_amp'], 0.0)
        * max(r['dominance_over_axis'], 0.0)
        * max(r['dominance_over_pair'], 0.0)
    )


def analyze(beta: float, ncell: int, pts_per_cell: int, modes: int, radius: float, skip_ground=True):
    mesh = build_supercell_mesh(ncell, pts_per_cell)
    A, M = rce.assemble(mesh, beta)
    vals, vecs = rce.solve_modes(A, M, modes)
    start = 1 if skip_ground else 0
    rows = [mode_metrics(mesh, vals, vecs, i, ncell, radius) for i in range(start, vecs.shape[1])]
    for r in rows:
        r['score'] = score_row(r)
        r['beta'] = beta
        r['ncell'] = ncell
        r['pts_per_cell'] = pts_per_cell
    best = max(rows, key=lambda r: r['score'])
    return mesh, vals, vecs, rows, best, rce.symmetry_norm(A)


def print_summary(beta, rows, best, max_rows=12):
    print(f'=== primitive cubic supercell ordered readout: beta={beta} ===')
    print('mode | lambda | mean|xyz| | uniform | best-q | q-amp | corr(x,y,z) | axis | pair | score')
    for r in rows[:max_rows]:
        qtxt = f"{r['best_q_label']}={r['best_q']}"
        ctxt = f"({r['corr_x']:+.2f},{r['corr_y']:+.2f},{r['corr_z']:+.2f})"
        print(f"{r['mode']:4d} | {r['lambda']:.6f} | {r['xyz_mean_abs']:.5f} | {r['xyz_uniformity']:.3f} | "
              f"{qtxt:16s} | {r['best_q_amp']:.3f} | {ctxt:18s} | {r['axis_mean_abs']:.3f} | {r['pair_mean_abs']:.3f} | {r['score']:.6f}")
    print('best nontrivial ordered xyz mode:')
    print(best)


def main():
    ap = argparse.ArgumentParser(description='Improved primitive cubic supercell readout with lattice-order diagnostics.')
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5,10')
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--radius', type=float, default=0.35)
    ap.add_argument('--include-ground', action='store_true')
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = [float(x) for x in args.betas.split(',') if x.strip()]
    allbest = []
    for beta in betas:
        mesh, vals, vecs, rows, best, sym = analyze(beta, args.ncell, args.pts_per_cell, args.modes, args.radius,
                                                    skip_ground=not args.include_ground)
        print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
        print(f'symmetry ||A-A^T||_F = {sym:.6e}')
        print_summary(beta, rows, best, max_rows=min(12, len(rows)))
        print()
        allbest.append(best)

    if args.csv:
        fields = ['beta','mode','lambda','xyz_mean_abs','xyz_std_abs','xyz_uniformity','xyz_same_sign',
                  'best_q','best_q_label','best_q_amp','corr_x','corr_y','corr_z','corr_mean',
                  'axis_mean_abs','pair_mean_abs','dominance_over_axis','dominance_over_pair','score','ncell','pts_per_cell']
        with open(args.csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in allbest:
                row = dict(r)
                row['best_q'] = str(row['best_q'])
                w.writerow({k: row[k] for k in fields})

if __name__ == '__main__':
    main()
