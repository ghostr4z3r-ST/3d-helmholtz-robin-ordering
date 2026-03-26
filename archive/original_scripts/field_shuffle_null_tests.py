import argparse, csv, importlib.util
from pathlib import Path
import numpy as np

spec = importlib.util.spec_from_file_location('ext', '/mnt/data/face_information_sidebranch_extended.py')
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)

FACES = ext.FACES
rng = np.random.default_rng(12345)


def collect_cube_modes(betas=(0.0,1.0,3.0,5.0), nx=17, ny=17, nz=17, modes=14, center_radius=0.22):
    items = []
    corners = ext.fis.rce.make_corners('cube')
    mesh = ext.fis.rce.build_hexahedral_mesh(nx, ny, nz, corners)
    for beta in betas:
        A, M = ext.fis.rce.assemble(mesh, beta)
        vals, vecs = ext.fis.rce.solve_modes(A, M, modes)
        sets = ext.fis.select_mode_sets(mesh, vals, vecs, beta, center_radius=center_radius)
        for label, idxs in sets.items():
            for mi in idxs:
                items.append({'mesh': mesh, 'beta': beta, 'label': label, 'mode_index': mi, 'mode': vecs[:, mi].copy()})
    return items


def features_for_mode(mesh, mode, face_radius=0.24, face_thickness=0.08, center_radius=0.22):
    face_sigs = {f: ext.fis.local_face_signature(mesh, mode, f, radius=face_radius, thickness=face_thickness) for f in FACES}
    onef = [{'label': None, 'feature': ext.improved_face_vector(face_sigs[f])} for f in FACES]
    sixf = ext.improved_sixface_vector(face_sigs)
    vol = ext.improved_volume_vector(mesh, mode, center_radius=center_radius)
    return onef, sixf, vol


def build_samples(items, shuffle=False):
    s1, s6, sv = [], [], []
    for it in items:
        mode = it['mode']
        if shuffle:
            mode = mode[rng.permutation(len(mode))]
        onef, sixf, vol = features_for_mode(it['mesh'], mode)
        for d in onef:
            d['label'] = it['label']
            s1.append(d)
        s6.append({'label': it['label'], 'feature': sixf})
        sv.append({'label': it['label'], 'feature': vol})
    return {'cube_raw_1F': s1, 'cube_raw_6F': s6, 'cube_volume': sv}


def accs(samples_dict):
    return {k: ext.nearest_centroid_cv(v) for k, v in samples_dict.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n_perm', type=int, default=200)
    ap.add_argument('--csv', type=str, default='/mnt/data/field_shuffle_null_tests.csv')
    args = ap.parse_args()

    items = collect_cube_modes()
    real = accs(build_samples(items, shuffle=False))
    perms = {k: [] for k in real}
    for _ in range(args.n_perm):
        sh = accs(build_samples(items, shuffle=True))
        for k,v in sh.items():
            perms[k].append(v)

    rows=[]
    for k, ar in real.items():
        arr = np.asarray(perms[k], float)
        rows.append({
            'task': k,
            'A_real': float(ar),
            'shuffle_mean': float(arr.mean()),
            'shuffle_std': float(arr.std(ddof=1)),
            'z_score': float((ar - arr.mean()) / (arr.std(ddof=1) + 1e-15)),
            'p_ge': float((np.sum(arr >= ar) + 1) / (len(arr) + 1)),
            'n_perm': int(args.n_perm),
        })
    out = Path(args.csv)
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(out)
    for r in rows:
        print(r)

if __name__ == '__main__':
    main()
