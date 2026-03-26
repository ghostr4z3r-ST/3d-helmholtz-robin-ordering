import argparse, importlib.util, csv, math
from pathlib import Path
import numpy as np

SRC = Path('/mnt/data/primitive_cubic_supercell_ordered.py')
spec = importlib.util.spec_from_file_location('pcs', SRC)
pcs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pcs)


def body_centers(ncell: int):
    return [(i + 0.5, j + 0.5, k + 0.5) for k in range(ncell) for j in range(ncell) for i in range(ncell)]


def interior_corners(ncell: int):
    if ncell < 2:
        return []
    return [(float(i), float(j), float(k))
            for k in range(1, ncell)
            for j in range(1, ncell)
            for i in range(1, ncell)]


def reshape_sites(values: np.ndarray, n: int):
    return values.reshape((n, n, n)) if n > 0 else np.zeros((0, 0, 0), dtype=float)


def site_metrics(mesh, modevec, sites, arr_n: int, radius_phys: float):
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
            'mean_abs': 0.0,
            'std_abs': 0.0,
            'uniformity': 0.0,
            'same_sign': 0.0,
            'best_q': (0, 0, 0),
            'best_q_label': 'none',
            'best_q_amp': 0.0,
            'corr_x': 0.0,
            'corr_y': 0.0,
            'corr_z': 0.0,
            'corr_mean': 0.0,
            'axis_mean_abs': 0.0,
            'pair_mean_abs': 0.0,
        }
    cube = reshape_sites(coeffs_xyz, arr_n)
    mean_abs = float(np.mean(np.abs(coeffs_xyz)))
    std_abs = float(np.std(np.abs(coeffs_xyz)))
    uniformity = float(max(0.0, 1.0 - std_abs / (mean_abs + 1e-15)))
    same_sign = float(abs(np.sum(coeffs_xyz)) / (np.sum(np.abs(coeffs_xyz)) + 1e-15))
    best_q, best_q_amp, _ = pcs.lattice_fourier_diagnostics(cube, arr_n)
    corr_x = pcs.neighbor_corr(cube, 0)
    corr_y = pcs.neighbor_corr(cube, 1)
    corr_z = pcs.neighbor_corr(cube, 2)
    return {
        'mean_abs': mean_abs,
        'std_abs': std_abs,
        'uniformity': uniformity,
        'same_sign': same_sign,
        'best_q': best_q,
        'best_q_label': pcs.classify_q(best_q, arr_n),
        'best_q_amp': best_q_amp,
        'corr_x': corr_x,
        'corr_y': corr_y,
        'corr_z': corr_z,
        'corr_mean': float((corr_x + corr_y + corr_z) / 3.0),
        'axis_mean_abs': float(np.mean(coeffs_axis)),
        'pair_mean_abs': float(np.mean(coeffs_pair)),
    }


def mode_metrics(mesh, vals, vecs, mode_index: int, ncell: int, radius: float):
    body = site_metrics(mesh, vecs[:, mode_index], body_centers(ncell), ncell, radius)
    corner = site_metrics(mesh, vecs[:, mode_index], interior_corners(ncell), max(ncell - 1, 0), radius)
    row = {
        'mode': mode_index + 1,
        'lambda': float(vals[mode_index]),
        'body_xyz_mean_abs': body['mean_abs'],
        'body_xyz_uniformity': body['uniformity'],
        'body_xyz_same_sign': body['same_sign'],
        'body_best_q': body['best_q'],
        'body_best_q_label': body['best_q_label'],
        'body_best_q_amp': body['best_q_amp'],
        'body_corr_x': body['corr_x'],
        'body_corr_y': body['corr_y'],
        'body_corr_z': body['corr_z'],
        'body_corr_mean': body['corr_mean'],
        'body_axis_mean_abs': body['axis_mean_abs'],
        'body_pair_mean_abs': body['pair_mean_abs'],
        'corner_xyz_mean_abs': corner['mean_abs'],
        'corner_xyz_uniformity': corner['uniformity'],
        'corner_xyz_same_sign': corner['same_sign'],
        'corner_best_q': corner['best_q'],
        'corner_best_q_label': corner['best_q_label'],
        'corner_best_q_amp': corner['best_q_amp'],
        'corner_corr_x': corner['corr_x'],
        'corner_corr_y': corner['corr_y'],
        'corner_corr_z': corner['corr_z'],
        'corner_corr_mean': corner['corr_mean'],
        'corner_axis_mean_abs': corner['axis_mean_abs'],
        'corner_pair_mean_abs': corner['pair_mean_abs'],
    }
    row['body_dom_axis'] = row['body_xyz_mean_abs'] / (row['body_axis_mean_abs'] + 1e-15)
    row['body_dom_pair'] = row['body_xyz_mean_abs'] / (row['body_pair_mean_abs'] + 1e-15)
    row['corner_dom_axis'] = row['corner_xyz_mean_abs'] / (row['corner_axis_mean_abs'] + 1e-15)
    row['corner_dom_pair'] = row['corner_xyz_mean_abs'] / (row['corner_pair_mean_abs'] + 1e-15)
    row['body_score'] = (
        row['body_xyz_mean_abs']
        * row['body_xyz_uniformity']
        * row['body_best_q_amp']
        * row['body_dom_axis']
        * row['body_dom_pair']
    )
    row['bcc_score'] = (
        math.sqrt(max(row['body_xyz_mean_abs'], 0.0) * max(row['corner_xyz_mean_abs'], 0.0))
        * math.sqrt(max(row['body_xyz_uniformity'], 0.0) * max(row['corner_xyz_uniformity'], 0.0))
        * math.sqrt(max(row['body_best_q_amp'], 0.0) * max(row['corner_best_q_amp'], 0.0))
    )
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
    best_body = max(rows, key=lambda r: r['body_score'])
    best_bcc = max(rows, key=lambda r: r['bcc_score'])
    sym = pcs.rce.symmetry_norm(A)
    return mesh, vals, vecs, rows, best_body, best_bcc, sym


def print_summary(beta, rows, best_body, best_bcc, max_rows=12):
    print(f'=== bcc-style readout on cubic supercell: beta={beta} ===')
    print('mode | lambda | body<|xyz|> | body-q | corner<|xyz|> | corner-q | body-score | bcc-score')
    for r in rows[:max_rows]:
        bq = f"{r['body_best_q_label']}={r['body_best_q']}"
        cq = f"{r['corner_best_q_label']}={r['corner_best_q']}"
        print(f"{r['mode']:4d} | {r['lambda']:.6f} | {r['body_xyz_mean_abs']:.5f} | {bq:16s} | {r['corner_xyz_mean_abs']:.5f} | {cq:16s} | {r['body_score']:.6f} | {r['bcc_score']:.6f}")
    print('best body-centered xyz mode:')
    print(best_body)
    print('best combined bcc-like body+corner mode:')
    print(best_bcc)


def main():
    ap = argparse.ArgumentParser(description='bcc-style diagnostic readout on cubic supercells.')
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
        mesh, vals, vecs, rows, best_body, best_bcc, sym = analyze(beta, args.ncell, args.pts_per_cell, args.modes, args.radius,
                                                                   skip_ground=not args.include_ground)
        print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
        print(f'symmetry ||A-A^T||_F = {sym:.6e}')
        print_summary(beta, rows, best_body, best_bcc, max_rows=min(len(rows), 12))
        print()
        row = dict(best_bcc)
        row['best_body_mode'] = best_body['mode']
        row['best_body_lambda'] = best_body['lambda']
        row['best_body_xyz_mean_abs'] = best_body['body_xyz_mean_abs']
        best_rows.append(row)
    if args.csv:
        fields = [
            'beta','ncell','pts_per_cell','mode','lambda',
            'body_xyz_mean_abs','body_xyz_uniformity','body_xyz_same_sign','body_best_q','body_best_q_label','body_best_q_amp',
            'corner_xyz_mean_abs','corner_xyz_uniformity','corner_xyz_same_sign','corner_best_q','corner_best_q_label','corner_best_q_amp',
            'body_score','bcc_score','best_body_mode','best_body_lambda','best_body_xyz_mean_abs'
        ]
        with open(args.csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in best_rows:
                rr = dict(r)
                rr['body_best_q'] = str(rr['body_best_q'])
                rr['corner_best_q'] = str(rr['corner_best_q'])
                w.writerow({k: rr.get(k, '') for k in fields})

if __name__ == '__main__':
    main()
