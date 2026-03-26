import argparse
import csv
import importlib.util
import math
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

SRC = Path('/mnt/data/robin_3d_center_emergence.py')
spec = importlib.util.spec_from_file_location('rce', SRC)
rce = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rce)

# ---------- local diagnostics ----------

def hadamard_basis_4():
    # quadrants in order (--,-+,+-,++) for two in-plane signs
    sigs = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]], dtype=int)
    names = ['const', 'a', 'b', 'ab']
    basis = []
    for name in names:
        if name == 'const':
            v = np.ones(4)
        elif name == 'a':
            v = sigs[:, 0]
        elif name == 'b':
            v = sigs[:, 1]
        elif name == 'ab':
            v = sigs[:, 0] * sigs[:, 1]
        v = v.astype(float)
        v /= np.linalg.norm(v)
        basis.append(v)
    return names, np.vstack(basis), sigs


def local_face_signature(mesh, mode, face: str, radius=0.22, thickness=0.10):
    lc = mesh.logical_coords
    names, basis, sigs = hadamard_basis_4()

    # choose face plane and in-plane coords
    if face == 'x0':
        d = lc[:, 0]
        in1, in2 = lc[:, 1] - 0.5, lc[:, 2] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    elif face == 'x1':
        d = 1.0 - lc[:, 0]
        in1, in2 = lc[:, 1] - 0.5, lc[:, 2] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    elif face == 'y0':
        d = lc[:, 1]
        in1, in2 = lc[:, 0] - 0.5, lc[:, 2] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    elif face == 'y1':
        d = 1.0 - lc[:, 1]
        in1, in2 = lc[:, 0] - 0.5, lc[:, 2] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    elif face == 'z0':
        d = lc[:, 2]
        in1, in2 = lc[:, 0] - 0.5, lc[:, 1] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    elif face == 'z1':
        d = 1.0 - lc[:, 2]
        in1, in2 = lc[:, 0] - 0.5, lc[:, 1] - 0.5
        mask = (d <= thickness) & (np.abs(in1) <= radius) & (np.abs(in2) <= radius)
    else:
        raise ValueError(face)

    eps = 1e-12
    mask &= (np.abs(in1) > eps) & (np.abs(in2) > eps)
    idxs = np.where(mask)[0]
    if len(idxs) == 0:
        raise RuntimeError(f'No points on face {face}; increase resolution/radius/thickness')

    vals4 = np.zeros(4, dtype=float)
    counts = np.zeros(4, dtype=int)
    for idx in idxs:
        s = np.array([1 if in1[idx] > 0 else -1, 1 if in2[idx] > 0 else -1], dtype=int)
        q = np.where((sigs == s).all(axis=1))[0][0]
        vals4[q] += mode[idx]
        counts[q] += 1
    if np.any(counts == 0):
        raise RuntimeError(f'Empty quadrants on face {face}: counts={counts.tolist()}')
    vals4 /= counts
    nrm = np.linalg.norm(vals4)
    vals4n = vals4 / nrm if nrm > 0 else vals4
    coeffs = basis @ vals4n
    order = np.argsort(-np.abs(coeffs))
    ranked = [(names[i], float(coeffs[i])) for i in order]
    return {
        'values4': vals4,
        'counts': counts,
        'coeffs': {names[i]: float(coeffs[i]) for i in range(4)},
        'ranked': ranked,
    }


def top_abs_coeff_name(coeffs: Dict[str, float]) -> str:
    return max(coeffs.items(), key=lambda kv: abs(kv[1]))[0]


# ---------- mode selection ----------

def center_signature(mesh, vec, radius=0.22):
    sig = rce.local_octant_signature(mesh, vec, (0.5, 0.5, 0.5), radius=radius)
    return sig['coeffs']


def select_mode_sets(mesh, vals, vecs, beta, center_radius=0.22):
    axis_candidates = []
    pair_candidates = []
    xyz_candidates = []
    for i in range(vecs.shape[1]):
        coeffs = center_signature(mesh, vecs[:, i], radius=center_radius)
        # skip ground-like const mode
        top = top_abs_coeff_name(coeffs)
        axis_score = max(abs(coeffs.get(k, 0.0)) for k in ['x', 'y', 'z'])
        pair_score = max(abs(coeffs.get(k, 0.0)) for k in ['xy', 'xz', 'yz'])
        xyz_score = abs(coeffs.get('xyz', 0.0))
        if top in {'x', 'y', 'z'}:
            axis_candidates.append((axis_score, i))
        if top in {'xy', 'xz', 'yz'}:
            pair_candidates.append((pair_score, i))
        if top == 'xyz':
            xyz_candidates.append((xyz_score, i))
    axis_candidates.sort(reverse=True)
    pair_candidates.sort(reverse=True)
    xyz_candidates.sort(reverse=True)
    axis_modes = sorted([i for _, i in axis_candidates[:3]])
    pair_modes = sorted([i for _, i in pair_candidates[:3]])
    xyz_mode = xyz_candidates[0][1] if xyz_candidates else None
    return {
        'axis': axis_modes,
        'pair': pair_modes,
        'full3d': [xyz_mode] if xyz_mode is not None else [],
    }


# ---------- features and simple classifiers ----------

def face_feature_vector(sig):
    coeffs = sig['coeffs']
    return np.array([coeffs['const'], coeffs['a'], coeffs['b'], coeffs['ab']], dtype=float)


def sixface_feature_vector(face_sigs: Dict[str, dict]):
    # histogram of top types + mean abs coeffs
    top_types = [top_abs_coeff_name(sig['coeffs']) for sig in face_sigs.values()]
    counts = Counter(top_types)
    hist = np.array([counts.get(k, 0) for k in ['const', 'a', 'b', 'ab']], dtype=float) / 6.0
    mean_abs = np.mean([
        [abs(sig['coeffs'][k]) for k in ['const', 'a', 'b', 'ab']]
        for sig in face_sigs.values()
    ], axis=0)
    return np.concatenate([hist, mean_abs])


def volume_feature_vector(mesh, vec, center_radius=0.22):
    coeffs = center_signature(mesh, vec, radius=center_radius)
    return np.array([coeffs[k] for k in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']], dtype=float)


def nearest_centroid_cv(samples: List[dict], feature_key: str) -> float:
    correct = 0
    for i, s in enumerate(samples):
        train = [t for j, t in enumerate(samples) if j != i]
        labels = sorted(set(t['label'] for t in train))
        centroids = {}
        for lab in labels:
            xs = np.vstack([t[feature_key] for t in train if t['label'] == lab])
            centroids[lab] = xs.mean(axis=0)
        x = s[feature_key]
        pred = min(labels, key=lambda lab: np.linalg.norm(x - centroids[lab]))
        if pred == s['label']:
            correct += 1
    return correct / len(samples) if samples else float('nan')


def family_from_face_types(face_types: List[str]) -> str:
    # coarse rule useful for interpretation only
    n_const = sum(t == 'const' for t in face_types)
    n_ab = sum(t == 'ab' for t in face_types)
    n_1d = 6 - n_const - n_ab
    return f"const={n_const},1d={n_1d},2d={n_ab}"


def run_first_test(betas=(0.0, 1.0, 5.0), nx=13, ny=13, nz=13, modes=12, face_radius=0.22, face_thickness=0.10, center_radius=0.22):
    corners = rce.make_corners('cube')
    mesh = rce.build_hexahedral_mesh(nx, ny, nz, corners)
    summary_rows = []
    face_samples = []
    sixface_samples = []
    volume_samples = []
    faces = ['x0', 'x1', 'y0', 'y1', 'z0', 'z1']

    for beta in betas:
        A, M = rce.assemble(mesh, beta)
        vals, vecs = rce.solve_modes(A, M, modes)
        sets = select_mode_sets(mesh, vals, vecs, beta, center_radius=center_radius)
        for label, idxs in sets.items():
            for mi in idxs:
                mode = vecs[:, mi]
                face_sigs = {f: local_face_signature(mesh, mode, f, radius=face_radius, thickness=face_thickness) for f in faces}
                ftypes = [top_abs_coeff_name(face_sigs[f]['coeffs']) for f in faces]
                summary_rows.append({
                    'beta': beta,
                    'label': label,
                    'mode': mi + 1,
                    'lambda': float(vals[mi]),
                    'face_pattern': '|'.join(ftypes),
                    'face_pattern_counts': family_from_face_types(ftypes),
                    'top_center_coeff': top_abs_coeff_name(center_signature(mesh, mode, radius=center_radius)),
                })
                for f in faces:
                    face_samples.append({
                        'label': label,
                        'feature': face_feature_vector(face_sigs[f]),
                    })
                sixface_samples.append({
                    'label': label,
                    'feature': sixface_feature_vector(face_sigs),
                })
                volume_samples.append({
                    'label': label,
                    'feature': volume_feature_vector(mesh, mode, center_radius=center_radius),
                })
    A1 = nearest_centroid_cv(face_samples, 'feature')
    A6 = nearest_centroid_cv(sixface_samples, 'feature')
    AV = nearest_centroid_cv(volume_samples, 'feature')
    D1 = 1.0 - (A1 / AV) if AV > 0 else float('nan')
    D6 = 1.0 - (A6 / AV) if AV > 0 else float('nan')
    return summary_rows, {'A_1F': A1, 'A_6F': A6, 'A_V': AV, 'D_1F': D1, 'D_6F': D6}


def main():
    ap = argparse.ArgumentParser(description='First test for face information deficit in cube Helmholtz-Robin modes.')
    ap.add_argument('--betas', type=str, default='0,1,5')
    ap.add_argument('--nx', type=int, default=13)
    ap.add_argument('--ny', type=int, default=13)
    ap.add_argument('--nz', type=int, default=13)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--csv', type=str, default='/mnt/data/face_information_first_test.csv')
    args = ap.parse_args()

    betas = tuple(float(x) for x in args.betas.split(','))
    rows, metrics = run_first_test(betas=betas, nx=args.nx, ny=args.ny, nz=args.nz, modes=args.modes)

    csv_path = Path(args.csv)
    with csv_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['beta', 'label', 'mode', 'lambda', 'top_center_coeff', 'face_pattern', 'face_pattern_counts'])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
        writer.writerow({})
        writer.writerow({'beta': 'METRICS', 'label': 'A_1F', 'mode': metrics['A_1F']})
        writer.writerow({'beta': 'METRICS', 'label': 'A_6F', 'mode': metrics['A_6F']})
        writer.writerow({'beta': 'METRICS', 'label': 'A_V', 'mode': metrics['A_V']})
        writer.writerow({'beta': 'METRICS', 'label': 'D_1F', 'mode': metrics['D_1F']})
        writer.writerow({'beta': 'METRICS', 'label': 'D_6F', 'mode': metrics['D_6F']})

    print('First test: face information deficit on cube Helmholtz-Robin modes')
    for r in rows:
        print(f"beta={r['beta']:>3} | {r['label']:>6} | mode={r['mode']:>2} | lambda={r['lambda']:.6f} | center={r['top_center_coeff']:<4} | faces={r['face_pattern_counts']} [{r['face_pattern']}]")
    print('\nMetrics:')
    for k, v in metrics.items():
        print(f'  {k} = {v:.6f}')
    print(f'CSV written to {csv_path}')

if __name__ == '__main__':
    main()
