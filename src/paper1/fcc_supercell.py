import argparse, csv, math
from pathlib import Path
import numpy as np

from paper1 import primitive_cubic_supercell as pcs


def cell_centers(ncell: int):
    sites = [(i + 0.5, j + 0.5, k + 0.5) for k in range(ncell) for j in range(ncell) for i in range(ncell)]
    shape = (ncell, ncell, ncell)  # z,y,x
    return sites, shape


def face_xy_sites(ncell: int):
    # interior xy faces between neighboring cells, located at integer z = 1..ncell-1
    sites = [(i + 0.5, j + 0.5, float(k)) for k in range(1, ncell) for j in range(ncell) for i in range(ncell)]
    shape = (max(ncell - 1, 0), ncell, ncell)
    return sites, shape


def face_xz_sites(ncell: int):
    sites = [(i + 0.5, float(j), k + 0.5) for k in range(ncell) for j in range(1, ncell) for i in range(ncell)]
    shape = (ncell, max(ncell - 1, 0), ncell)
    return sites, shape


def face_yz_sites(ncell: int):
    sites = [(float(i), j + 0.5, k + 0.5) for k in range(ncell) for j in range(ncell) for i in range(1, ncell)]
    shape = (ncell, ncell, max(ncell - 1, 0))
    return sites, shape


def lattice_fourier_rect(arr: np.ndarray):
    shape = arr.shape
    total_abs = float(np.sum(np.abs(arr))) + 1e-15
    best_amp = -1.0
    best_q = (0, 0, 0)
    q_amps = {}
    nz, ny, nx = shape
    for mz in range(nz):
        for my in range(ny):
            for mx in range(nx):
                s = 0.0 + 0.0j
                for k in range(nz):
                    for j in range(ny):
                        for i in range(nx):
                            phase = -2j * math.pi * ((mx * i / max(nx,1)) + (my * j / max(ny,1)) + (mz * k / max(nz,1)))
                            s += arr[k, j, i] * np.exp(phase)
                amp = abs(s) / total_abs
                q = (mx, my, mz)
                q_amps[q] = float(amp)
                if amp > best_amp:
                    best_amp = float(amp)
                    best_q = q
    return best_q, best_amp, q_amps


def neighbor_corr_rect(arr: np.ndarray, axis: int):
    total = 0.0
    denom = 0.0
    nz, ny, nx = arr.shape
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                a = arr[k, j, i]
                if axis == 0 and nx > 1:
                    b = arr[k, j, (i + 1) % nx]
                elif axis == 1 and ny > 1:
                    b = arr[k, (j + 1) % ny, i]
                elif axis == 2 and nz > 1:
                    b = arr[(k + 1) % nz, j, i]
                else:
                    continue
                total += a * b
                denom += abs(a * b)
    return 0.0 if denom < 1e-15 else float(total / denom)


def classify_q_rect(q, shape):
    labels = ['X', 'Y', 'Z']
    comps = [q[0], q[1], q[2]]
    dims = [shape[2], shape[1], shape[0]]
    parts = []
    for comp, dim, lab in zip(comps, dims, labels):
        if dim > 0 and comp % dim != 0:
            parts.append(f'{lab}{comp}')
    return '+'.join(parts) if parts else 'const'


def family_metrics(mesh, modevec, site_shape, radius_phys: float):
    sites, shape = site_shape
    coeffs_xyz = []
    coeffs_axis = []
    coeffs_pair = []
    for c in sites:
        sig = pcs.local_sig_physical(mesh, modevec, c, radius_phys)
        coeffs = sig['coeffs']
        coeffs_xyz.append(coeffs['xyz'])
        coeffs_axis.append(max(abs(coeffs['x']), abs(coeffs['y']), abs(coeffs['z'])))
        coeffs_pair.append(max(abs(coeffs['xy']), abs(coeffs['xz']), abs(coeffs['yz'])))
    coeffs_xyz = np.asarray(coeffs_xyz, dtype=float)
    coeffs_axis = np.asarray(coeffs_axis, dtype=float)
    coeffs_pair = np.asarray(coeffs_pair, dtype=float)
    if coeffs_xyz.size == 0:
        return {
            'xyz_mean_abs': 0.0, 'xyz_uniformity': 0.0, 'xyz_same_sign': 0.0,
            'best_q': (0,0,0), 'best_q_label': 'none', 'best_q_amp': 0.0,
            'corr_x': 0.0, 'corr_y': 0.0, 'corr_z': 0.0, 'corr_mean': 0.0,
            'axis_mean_abs': 0.0, 'pair_mean_abs': 0.0, 'dom_axis': 0.0, 'dom_pair': 0.0,
        }
    arr = coeffs_xyz.reshape(shape)
    mean_abs = float(np.mean(np.abs(coeffs_xyz)))
    std_abs = float(np.std(np.abs(coeffs_xyz)))
    uniformity = float(max(0.0, 1.0 - std_abs / (mean_abs + 1e-15)))
    same_sign = float(abs(np.sum(coeffs_xyz)) / (np.sum(np.abs(coeffs_xyz)) + 1e-15))
    best_q, best_q_amp, _ = lattice_fourier_rect(arr)
    corr_x = neighbor_corr_rect(arr, 0)
    corr_y = neighbor_corr_rect(arr, 1)
    corr_z = neighbor_corr_rect(arr, 2)
    return {
        'xyz_mean_abs': mean_abs,
        'xyz_uniformity': uniformity,
        'xyz_same_sign': same_sign,
        'best_q': best_q,
        'best_q_label': classify_q_rect(best_q, shape),
        'best_q_amp': float(best_q_amp),
        'corr_x': corr_x,
        'corr_y': corr_y,
        'corr_z': corr_z,
        'corr_mean': float((corr_x + corr_y + corr_z) / 3.0),
        'axis_mean_abs': float(np.mean(coeffs_axis)),
        'pair_mean_abs': float(np.mean(coeffs_pair)),
        'dom_axis': float(mean_abs / (float(np.mean(coeffs_axis)) + 1e-15)),
        'dom_pair': float(mean_abs / (float(np.mean(coeffs_pair)) + 1e-15)),
    }


def mode_metrics(mesh, vals, vecs, mode_index: int, ncell: int, radius_phys: float):
    v = vecs[:, mode_index]
    center = family_metrics(mesh, v, cell_centers(ncell), radius_phys)
    xy = family_metrics(mesh, v, face_xy_sites(ncell), radius_phys)
    xz = family_metrics(mesh, v, face_xz_sites(ncell), radius_phys)
    yz = family_metrics(mesh, v, face_yz_sites(ncell), radius_phys)

    row = {
        'mode': mode_index + 1,
        'lambda': float(vals[mode_index]),
        'center_xyz_mean_abs': center['xyz_mean_abs'],
        'center_xyz_uniformity': center['xyz_uniformity'],
        'center_best_q': center['best_q'],
        'center_best_q_label': center['best_q_label'],
        'center_best_q_amp': center['best_q_amp'],
        'xy_xyz_mean_abs': xy['xyz_mean_abs'],
        'xy_xyz_uniformity': xy['xyz_uniformity'],
        'xy_best_q': xy['best_q'],
        'xy_best_q_label': xy['best_q_label'],
        'xy_best_q_amp': xy['best_q_amp'],
        'xz_xyz_mean_abs': xz['xyz_mean_abs'],
        'xz_xyz_uniformity': xz['xyz_uniformity'],
        'xz_best_q': xz['best_q'],
        'xz_best_q_label': xz['best_q_label'],
        'xz_best_q_amp': xz['best_q_amp'],
        'yz_xyz_mean_abs': yz['xyz_mean_abs'],
        'yz_xyz_uniformity': yz['xyz_uniformity'],
        'yz_best_q': yz['best_q'],
        'yz_best_q_label': yz['best_q_label'],
        'yz_best_q_amp': yz['best_q_amp'],
        'xy_dom_axis': xy['dom_axis'], 'xy_dom_pair': xy['dom_pair'],
        'xz_dom_axis': xz['dom_axis'], 'xz_dom_pair': xz['dom_pair'],
        'yz_dom_axis': yz['dom_axis'], 'yz_dom_pair': yz['dom_pair'],
    }
    face_mean = (xy['xyz_mean_abs'] * xz['xyz_mean_abs'] * yz['xyz_mean_abs']) ** (1/3)
    face_uniform = (xy['xyz_uniformity'] * xz['xyz_uniformity'] * yz['xyz_uniformity']) ** (1/3)
    face_qamp = (xy['best_q_amp'] * xz['best_q_amp'] * yz['best_q_amp']) ** (1/3)
    face_dom_axis = (xy['dom_axis'] * xz['dom_axis'] * yz['dom_axis']) ** (1/3)
    face_dom_pair = (xy['dom_pair'] * xz['dom_pair'] * yz['dom_pair']) ** (1/3)
    row['fcc_face_mean'] = face_mean
    row['fcc_face_uniformity'] = face_uniform
    row['fcc_face_qamp'] = face_qamp
    row['fcc_face_dom_axis'] = face_dom_axis
    row['fcc_face_dom_pair'] = face_dom_pair
    row['fcc_score'] = face_mean * face_uniform * face_qamp * face_dom_axis * face_dom_pair
    row['center_score'] = center['xyz_mean_abs'] * center['xyz_uniformity'] * center['best_q_amp'] * center['dom_axis'] * center['dom_pair']
    return row


def analyze(beta: float, ncell: int, pts_per_cell: int, modes: int, radius: float, skip_ground=True):
    mesh = pcs.build_supercell_mesh(ncell, pts_per_cell)
    A, M = pcs.rce.assemble(mesh, beta)
    vals, vecs = pcs.rce.solve_modes(A, M, modes)
    start = 1 if skip_ground else 0
    rows = [mode_metrics(mesh, vals, vecs, i, ncell, radius) for i in range(start, vecs.shape[1])]
    for r in rows:
        r['beta'] = beta
        r['ncell'] = ncell
        r['pts_per_cell'] = pts_per_cell
    best_center = max(rows, key=lambda r: r['center_score'])
    best_fcc = max(rows, key=lambda r: r['fcc_score'])
    sym = pcs.rce.symmetry_norm(A)
    return mesh, vals, vecs, rows, best_center, best_fcc, sym


def print_summary(beta, rows, best_center, best_fcc, max_rows=12):
    print(f'=== fcc-style readout on cubic supercell: beta={beta} ===')
    print('mode | lambda | center<|xyz|> | xy<|xyz|> | xz<|xyz|> | yz<|xyz|> | center-score | fcc-score')
    for r in rows[:max_rows]:
        print(f"{r['mode']:4d} | {r['lambda']:.6f} | {r['center_xyz_mean_abs']:.5f} | {r['xy_xyz_mean_abs']:.5f} | {r['xz_xyz_mean_abs']:.5f} | {r['yz_xyz_mean_abs']:.5f} | {r['center_score']:.6f} | {r['fcc_score']:.6f}")
    print('best center-dominated xyz mode:')
    print(best_center)
    print('best combined fcc-like face-centered xyz mode:')
    print(best_fcc)


def main():
    ap = argparse.ArgumentParser(description='fcc-style diagnostic readout on cubic supercells.')
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5,10')
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--radius', type=float, default=0.28)
    ap.add_argument('--include-ground', action='store_true')
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = [float(x) for x in args.betas.split(',') if x.strip()]
    best_rows = []
    for beta in betas:
        mesh, vals, vecs, rows, best_center, best_fcc, sym = analyze(beta, args.ncell, args.pts_per_cell, args.modes, args.radius,
                                                                     skip_ground=not args.include_ground)
        print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
        print(f'symmetry ||A-A^T||_F = {sym:.6e}')
        print_summary(beta, rows, best_center, best_fcc, max_rows=min(len(rows), 12))
        print()
        row = dict(best_fcc)
        row['best_center_mode'] = best_center['mode']
        row['best_center_lambda'] = best_center['lambda']
        row['best_center_xyz_mean_abs'] = best_center['center_xyz_mean_abs']
        best_rows.append(row)
    if args.csv:
        out_path = Path(args.csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fields = [
            'beta','ncell','pts_per_cell','mode','lambda',
            'center_xyz_mean_abs','center_xyz_uniformity','center_best_q','center_best_q_label','center_best_q_amp',
            'xy_xyz_mean_abs','xy_xyz_uniformity','xy_best_q','xy_best_q_label','xy_best_q_amp',
            'xz_xyz_mean_abs','xz_xyz_uniformity','xz_best_q','xz_best_q_label','xz_best_q_amp',
            'yz_xyz_mean_abs','yz_xyz_uniformity','yz_best_q','yz_best_q_label','yz_best_q_amp',
            'fcc_face_mean','fcc_face_uniformity','fcc_face_qamp','fcc_face_dom_axis','fcc_face_dom_pair','fcc_score',
            'best_center_mode','best_center_lambda','best_center_xyz_mean_abs'
        ]
        with out_path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in best_rows:
                rr = dict(r)
                rr['center_best_q'] = str(rr['center_best_q'])
                rr['xy_best_q'] = str(rr['xy_best_q'])
                rr['xz_best_q'] = str(rr['xz_best_q'])
                rr['yz_best_q'] = str(rr['yz_best_q'])
                w.writerow({k: rr.get(k, '') for k in fields})

if __name__ == '__main__':
    main()
