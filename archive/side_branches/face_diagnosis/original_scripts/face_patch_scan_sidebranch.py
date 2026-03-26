import argparse
import csv
import importlib.util
from pathlib import Path
from collections import defaultdict
import numpy as np

# Imports from prior work
spec = importlib.util.spec_from_file_location('fis', '/mnt/data/face_information_sidebranch.py')
fis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fis)

spec2 = importlib.util.spec_from_file_location('fise', '/mnt/data/face_information_sidebranch_extended.py')
fise = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(fise)

FACES = ['x0','x1','y0','y1','z0','z1']


def face_grid(mesh, mode, face):
    lc = mesh.logical_coords
    # identify face nodes
    if face == 'x0':
        mask = np.isclose(lc[:,0], 0.0)
        a, b = lc[mask,1], lc[mask,2]
    elif face == 'x1':
        mask = np.isclose(lc[:,0], 1.0)
        a, b = lc[mask,1], lc[mask,2]
    elif face == 'y0':
        mask = np.isclose(lc[:,1], 0.0)
        a, b = lc[mask,0], lc[mask,2]
    elif face == 'y1':
        mask = np.isclose(lc[:,1], 1.0)
        a, b = lc[mask,0], lc[mask,2]
    elif face == 'z0':
        mask = np.isclose(lc[:,2], 0.0)
        a, b = lc[mask,0], lc[mask,1]
    elif face == 'z1':
        mask = np.isclose(lc[:,2], 1.0)
        a, b = lc[mask,0], lc[mask,1]
    else:
        raise ValueError(face)
    vals = mode[mask]
    au = np.unique(np.round(a, 12))
    bu = np.unique(np.round(b, 12))
    ia = {float(v): i for i, v in enumerate(au)}
    ib = {float(v): i for i, v in enumerate(bu)}
    G = np.zeros((len(au), len(bu)), dtype=float)
    for aa, bb, vv in zip(a, b, vals):
        G[ia[float(round(aa,12))], ib[float(round(bb,12))]] = vv
    return G


def patch_features(P):
    P = np.asarray(P, dtype=float)
    mean = float(P.mean())
    absmean = float(np.mean(np.abs(P)))
    std = float(P.std())
    sign_bal = float(abs(np.mean(np.sign(P + 1e-15))))
    # quadrant averages
    m, n = P.shape
    mx, my = m // 2, n // 2
    # allow odd patch sizes by including center in upper/right blocks
    quads = np.array([
        P[:mx, :my].mean(),
        P[:mx, my:].mean(),
        P[mx:, :my].mean(),
        P[mx:, my:].mean(),
    ], dtype=float)
    qn = quads / (np.linalg.norm(quads) + 1e-15)
    basis = {
        'const': np.array([1,1,1,1], dtype=float) / 2.0,
        'a': np.array([-1,-1,1,1], dtype=float) / 2.0,
        'b': np.array([-1,1,-1,1], dtype=float) / 2.0,
        'ab': np.array([1,-1,-1,1], dtype=float) / 2.0,
    }
    coeffs = {k: float(v @ qn) for k, v in basis.items()}
    # simple gradient/alternation proxies
    dx = float(np.mean(np.abs(np.diff(P, axis=0)))) if m > 1 else 0.0
    dy = float(np.mean(np.abs(np.diff(P, axis=1)))) if n > 1 else 0.0
    checker = np.fromfunction(lambda i, j: (-1.0) ** (i + j), P.shape)
    stripe_x = np.fromfunction(lambda i, j: (-1.0) ** i, P.shape)
    stripe_y = np.fromfunction(lambda i, j: (-1.0) ** j, P.shape)
    c_checker = float(np.sum(P * checker) / (np.linalg.norm(P) * np.linalg.norm(checker) + 1e-15))
    c_x = float(np.sum(P * stripe_x) / (np.linalg.norm(P) * np.linalg.norm(stripe_x) + 1e-15))
    c_y = float(np.sum(P * stripe_y) / (np.linalg.norm(P) * np.linalg.norm(stripe_y) + 1e-15))
    return np.array([
        mean, absmean, std, sign_bal,
        coeffs['const'], coeffs['a'], coeffs['b'], coeffs['ab'],
        abs(coeffs['const']), abs(coeffs['a']), abs(coeffs['b']), abs(coeffs['ab']),
        dx, dy, c_x, c_y, c_checker
    ], dtype=float)


def scan_face_feature(G, patch_sizes=(3,5)):
    all_feat = []
    for p in patch_sizes:
        if G.shape[0] < p or G.shape[1] < p:
            continue
        for i in range(G.shape[0] - p + 1):
            for j in range(G.shape[1] - p + 1):
                all_feat.append(patch_features(G[i:i+p, j:j+p]))
    arr = np.vstack(all_feat)
    # summarize patch statistics
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    q25 = np.quantile(arr, 0.25, axis=0)
    q75 = np.quantile(arr, 0.75, axis=0)
    return np.concatenate([mean, std, q25, q75])


def scan_sixface_feature(face_grids, patch_sizes=(3,5)):
    per_face = np.vstack([scan_face_feature(face_grids[f], patch_sizes=patch_sizes) for f in FACES])
    mean = per_face.mean(axis=0)
    std = per_face.std(axis=0)
    # opposite-face cosine similarities
    pairs = [('x0','x1'), ('y0','y1'), ('z0','z1')]
    cors = []
    feats = {f: scan_face_feature(face_grids[f], patch_sizes=patch_sizes) for f in FACES}
    for a, b in pairs:
        va, vb = feats[a], feats[b]
        cors.append(float(va @ vb / (np.linalg.norm(va) * np.linalg.norm(vb) + 1e-15)))
    return np.concatenate([mean, std, np.array(cors)])


def nearest_centroid_cv(samples, feature_key='feature'):
    if not samples:
        return float('nan')
    correct = 0
    for i, s in enumerate(samples):
        train = [t for j,t in enumerate(samples) if j != i]
        labels = sorted(set(t['label'] for t in train))
        centroids = {lab: np.vstack([t[feature_key] for t in train if t['label']==lab]).mean(axis=0) for lab in labels}
        x = s[feature_key]
        pred = min(labels, key=lambda lab: np.linalg.norm(x - centroids[lab]))
        correct += int(pred == s['label'])
    return correct / len(samples)


def run_test(betas=(0.0,1.0,5.0), nx=17, ny=17, nz=17, modes=14, center_radius=0.22, face_radius=0.24, face_thickness=0.08, patch_sizes=(3,5)):
    corners = fis.rce.make_corners('cube')
    mesh = fis.rce.build_hexahedral_mesh(nx, ny, nz, corners)
    rows = []
    raw_1f, raw_6f, scan_1f, scan_6f, vols = [], [], [], [], []
    for beta in betas:
        A, M = fis.rce.assemble(mesh, beta)
        vals, vecs = fis.rce.solve_modes(A, M, modes)
        sets = fis.select_mode_sets(mesh, vals, vecs, beta, center_radius=center_radius)
        for label, idxs in sets.items():
            for mi in idxs:
                mode = vecs[:, mi]
                face_sigs = {f: fis.local_face_signature(mesh, mode, f, radius=face_radius, thickness=face_thickness) for f in FACES}
                face_grids = {f: face_grid(mesh, mode, f) for f in FACES}
                for f in FACES:
                    raw_1f.append({'label': label, 'feature': fise.improved_face_vector(face_sigs[f])})
                    scan_1f.append({'label': label, 'feature': scan_face_feature(face_grids[f], patch_sizes=patch_sizes)})
                raw_6f.append({'label': label, 'feature': fise.improved_sixface_vector(face_sigs)})
                scan_6f.append({'label': label, 'feature': scan_sixface_feature(face_grids, patch_sizes=patch_sizes)})
                vols.append({'label': label, 'feature': fise.improved_volume_vector(mesh, mode, center_radius=center_radius)})
                rows.append({
                    'beta': beta,
                    'label': label,
                    'mode': mi + 1,
                    'lambda': float(vals[mi]),
                    'center_top': fis.top_abs_coeff_name(fis.center_signature(mesh, mode, radius=center_radius)),
                    'raw_face_pattern': '|'.join(fis.top_abs_coeff_name(face_sigs[f]['coeffs']) for f in FACES),
                    'grid_shape': f"{face_grids['x0'].shape[0]}x{face_grids['x0'].shape[1]}"
                })
    metrics = {}
    metrics['A_1F'] = nearest_centroid_cv(raw_1f)
    metrics['A_6F'] = nearest_centroid_cv(raw_6f)
    metrics['A_1F_scan'] = nearest_centroid_cv(scan_1f)
    metrics['A_6F_scan'] = nearest_centroid_cv(scan_6f)
    metrics['A_V'] = nearest_centroid_cv(vols)
    for k in ['A_1F', 'A_6F', 'A_1F_scan', 'A_6F_scan']:
        metrics['D_' + k[2:]] = 1.0 - metrics[k] / metrics['A_V'] if metrics['A_V'] > 0 else float('nan')
    metrics['G_1F'] = metrics['A_1F_scan'] - metrics['A_1F']
    metrics['G_6F'] = metrics['A_6F_scan'] - metrics['A_6F']
    return rows, metrics


def main():
    ap = argparse.ArgumentParser(description='Statistical face-scan sidebranch test on cube Helmholtz-Robin modes.')
    ap.add_argument('--betas', type=str, default='0,1,5')
    ap.add_argument('--nx', type=int, default=17)
    ap.add_argument('--ny', type=int, default=17)
    ap.add_argument('--nz', type=int, default=17)
    ap.add_argument('--modes', type=int, default=14)
    ap.add_argument('--patch-sizes', type=str, default='3,5')
    ap.add_argument('--csv', type=str, default='/mnt/data/face_patch_scan_first_test.csv')
    args = ap.parse_args()

    betas = tuple(float(x) for x in args.betas.split(','))
    patch_sizes = tuple(int(x) for x in args.patch_sizes.split(','))
    rows, metrics = run_test(betas=betas, nx=args.nx, ny=args.ny, nz=args.nz, modes=args.modes, patch_sizes=patch_sizes)

    out = Path(args.csv)
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['section','beta','label','mode','lambda','center_top','raw_face_pattern','grid_shape','metric','value'])
        writer.writeheader()
        for r in rows:
            writer.writerow({'section':'rows', **r, 'metric':'', 'value':''})
        for k,v in metrics.items():
            writer.writerow({'section':'metrics','beta':'','label':'','mode':'','lambda':'','center_top':'','raw_face_pattern':'','grid_shape':'','metric':k,'value':v})

    print('Statistical face-scan sidebranch test')
    for k,v in metrics.items():
        print(f'  {k} = {v}')
    print(f'CSV written to {out}')


if __name__ == '__main__':
    main()
