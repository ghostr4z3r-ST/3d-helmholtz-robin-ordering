import argparse, importlib.util, csv, math
from pathlib import Path
import numpy as np
import pandas as pd


def load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

pcs = load_module('pcs', '/mnt/data/primitive_cubic_supercell_ordered.py')
fcc = load_module('fcc', '/mnt/data/fcc_supercell_readout.py')
bcc = load_module('bcc', '/mnt/data/bcc_supercell_readout.py')


def structure_factor(positions: np.ndarray, amps: np.ndarray, ncell: int):
    """Normalized structure factor on integer reciprocal grid 0..ncell-1.
    S(h,k,l) = |sum_j a_j exp(-2pi i (h x + k y + l z)/ncell)|^2 / (sum |a_j|)^2
    """
    norm = float(np.sum(np.abs(amps))) + 1e-15
    grid = {}
    best_q = (0, 0, 0)
    best_s = -1.0
    total = 0.0
    count = 0
    for h in range(ncell):
        for k in range(ncell):
            for l in range(ncell):
                phase = np.exp(-2j * math.pi * (h * positions[:, 0] + k * positions[:, 1] + l * positions[:, 2]) / ncell)
                s = abs(np.sum(amps * phase)) ** 2 / (norm ** 2)
                q = (h, k, l)
                grid[q] = float(s)
                total += float(s)
                count += 1
                if q != (0, 0, 0) and s > best_s:
                    best_s = float(s)
                    best_q = q
    mean_s = total / max(count, 1)
    contrast = (best_s - mean_s) / (best_s + mean_s + 1e-15)
    q0 = grid[(0, 0, 0)]
    ratio = best_s / (q0 + 1e-15)
    return {
        'best_q': best_q,
        'best_q_label': pcs.classify_q(best_q, ncell),
        'best_S': best_s,
        'mean_S': mean_s,
        'contrast': contrast,
        'S0': q0,
        'peak_to_q0': ratio,
    }


def xyz_coeffs_at_sites(mesh, modevec, sites, radius):
    vals = []
    for site in sites:
        sig = pcs.local_sig_physical(mesh, modevec, site, radius)
        vals.append(sig['coeffs']['xyz'])
    return np.asarray(vals, dtype=float)


def primitive_family(mesh, vecs, mode_index, ncell, radius):
    sites = pcs.subcell_centers(ncell)
    amps = xyz_coeffs_at_sites(mesh, vecs[:, mode_index], sites, radius)
    sf = structure_factor(np.asarray(sites, dtype=float), amps, ncell)
    return {
        'family': 'primitive',
        'subfamily': 'centers',
        'site_count': len(sites),
        'xyz_mean_abs': float(np.mean(np.abs(amps))),
        'uniformity': float(max(0.0, 1.0 - np.std(np.abs(amps)) / (np.mean(np.abs(amps)) + 1e-15))),
        'same_sign': float(abs(np.sum(amps)) / (np.sum(np.abs(amps)) + 1e-15)),
        **sf,
    }


def fcc_family(mesh, vecs, mode_index, ncell, radius):
    xy_sites, _ = fcc.face_xy_sites(ncell)
    xz_sites, _ = fcc.face_xz_sites(ncell)
    yz_sites, _ = fcc.face_yz_sites(ncell)
    positions = np.asarray(xy_sites + xz_sites + yz_sites, dtype=float)
    amps = np.concatenate([
        xyz_coeffs_at_sites(mesh, vecs[:, mode_index], xy_sites, radius),
        xyz_coeffs_at_sites(mesh, vecs[:, mode_index], xz_sites, radius),
        xyz_coeffs_at_sites(mesh, vecs[:, mode_index], yz_sites, radius),
    ])
    sf = structure_factor(positions, amps, ncell)
    return {
        'family': 'fcc',
        'subfamily': 'faces_combined',
        'site_count': len(positions),
        'xyz_mean_abs': float(np.mean(np.abs(amps))),
        'uniformity': float(max(0.0, 1.0 - np.std(np.abs(amps)) / (np.mean(np.abs(amps)) + 1e-15))),
        'same_sign': float(abs(np.sum(amps)) / (np.sum(np.abs(amps)) + 1e-15)),
        **sf,
    }


def bcc_family(mesh, vecs, mode_index, ncell, radius):
    body_sites = bcc.body_centers(ncell)
    corner_sites = bcc.interior_corners(ncell)
    body_amps = xyz_coeffs_at_sites(mesh, vecs[:, mode_index], body_sites, radius)
    corner_amps = xyz_coeffs_at_sites(mesh, vecs[:, mode_index], corner_sites, radius) if corner_sites else np.array([], dtype=float)
    positions = np.asarray(body_sites + corner_sites, dtype=float)
    amps = np.concatenate([body_amps, corner_amps]) if corner_sites else body_amps.copy()
    sf = structure_factor(positions, amps, ncell)
    row = {
        'family': 'bcc',
        'subfamily': 'body_plus_corners',
        'site_count': len(positions),
        'xyz_mean_abs': float(np.mean(np.abs(amps))),
        'uniformity': float(max(0.0, 1.0 - np.std(np.abs(amps)) / (np.mean(np.abs(amps)) + 1e-15))),
        'same_sign': float(abs(np.sum(amps)) / (np.sum(np.abs(amps)) + 1e-15)),
        'body_mean_abs': float(np.mean(np.abs(body_amps))) if body_amps.size else 0.0,
        'corner_mean_abs': float(np.mean(np.abs(corner_amps))) if corner_amps.size else 0.0,
        **sf,
    }
    row['P_bc'] = (row['body_mean_abs'] - row['corner_mean_abs']) / (row['body_mean_abs'] + row['corner_mean_abs'] + 1e-15)
    return row


def analyze_beta(beta, ncell=3, pts_per_cell=6, modes=20, radius=0.28):
    mesh = pcs.build_supercell_mesh(ncell, pts_per_cell)
    A, M = pcs.rce.assemble(mesh, beta)
    vals, vecs = pcs.rce.solve_modes(A, M, modes)

    # get best family modes from existing diagnostics
    _, _, _, prow_rows, pbest, _ = pcs.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
    _, _, _, frows, best_center, best_fcc, _ = fcc.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
    _, _, _, brows, best_body, best_bcc, _ = bcc.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)

    families = [
        ('primitive', pbest['mode'] - 1, pbest['lambda']),
        ('fcc', best_fcc['mode'] - 1, best_fcc['lambda']),
        ('bcc', best_bcc['mode'] - 1, best_bcc['lambda']),
    ]

    rows = []
    for fam, idx, lam in families:
        if fam == 'primitive':
            row = primitive_family(mesh, vecs, idx, ncell, radius)
        elif fam == 'fcc':
            row = fcc_family(mesh, vecs, idx, ncell, radius)
        else:
            row = bcc_family(mesh, vecs, idx, ncell, radius)
        row['beta'] = beta
        row['mode'] = idx + 1
        row['lambda'] = float(lam)
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description='Structure-factor proxy for primitive cubic, bcc, and fcc readouts.')
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5,10')
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--radius', type=float, default=0.28)
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = [float(x) for x in args.betas.split(',') if x.strip()]
    all_rows = []
    for beta in betas:
        all_rows.extend(analyze_beta(beta, args.ncell, args.pts_per_cell, args.modes, args.radius))
    df = pd.DataFrame(all_rows)
    with pd.option_context('display.max_columns', None, 'display.width', 220):
        print(df[['beta','family','mode','lambda','xyz_mean_abs','uniformity','same_sign','best_q_label','best_q','best_S','contrast','peak_to_q0'] + ([ 'body_mean_abs','corner_mean_abs','P_bc'] if 'P_bc' in df.columns else [])])

    if args.csv:
        out = df.copy()
        out['best_q'] = out['best_q'].astype(str)
        out.to_csv(args.csv, index=False)


if __name__ == '__main__':
    main()
