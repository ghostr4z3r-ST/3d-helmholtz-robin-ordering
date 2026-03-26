
import argparse, csv, importlib.util
from pathlib import Path
import numpy as np

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

fis = load('fis', '/mnt/data/face_information_sidebranch.py')
fise = load('fise', '/mnt/data/face_information_sidebranch_extended.py')
fps = load('fps', '/mnt/data/face_patch_scan_sidebranch.py')
fpss = load('fpss', '/mnt/data/face_patch_scan_supercell.py')
pcs = load('pcs', '/mnt/data/primitive_cubic_supercell_ordered.py')
fccr = load('fccr', '/mnt/data/fcc_supercell_readout.py')

FACES = ['x0','x1','y0','y1','z0','z1']

def nearest_centroid_accuracy(features, labels):
    features = [np.asarray(x, dtype=float) for x in features]
    labels = list(labels)
    n = len(labels)
    if n < 2:
        return float('nan')
    correct = 0
    for i in range(n):
        train_idx = [j for j in range(n) if j != i]
        labs = sorted(set(labels[j] for j in train_idx))
        centroids = {}
        for lab in labs:
            xs = [features[j] for j in train_idx if labels[j] == lab]
            centroids[lab] = np.vstack(xs).mean(axis=0)
        x = features[i]
        pred = min(labs, key=lambda lab: np.linalg.norm(x - centroids[lab]))
        correct += int(pred == labels[i])
    return correct / n

def permutation_test(features, labels, n_perm=300, seed=0):
    rng = np.random.default_rng(seed)
    real = nearest_centroid_accuracy(features, labels)
    perm = np.empty(n_perm, dtype=float)
    labels = np.array(labels, dtype=object)
    for k in range(n_perm):
        shuf = labels.copy()
        rng.shuffle(shuf)
        perm[k] = nearest_centroid_accuracy(features, shuf.tolist())
    p_ge = (1 + np.sum(perm >= real)) / (n_perm + 1)
    z = (real - perm.mean()) / (perm.std(ddof=1) + 1e-15)
    return {
        'A_real': float(real),
        'A_perm_mean': float(perm.mean()),
        'A_perm_std': float(perm.std(ddof=1)),
        'A_perm_q95': float(np.quantile(perm, 0.95)),
        'A_perm_q99': float(np.quantile(perm, 0.99)),
        'p_ge': float(p_ge),
        'z_score': float(z),
        'n_samples': int(len(labels)),
        'n_classes': int(len(set(labels)))
    }

def build_cube_datasets(betas=(0.0,1.0,3.0,5.0), nx=17, ny=17, nz=17, modes=14, face_radius=0.24, face_thickness=0.08, center_radius=0.22):
    corners = fis.rce.make_corners('cube')
    mesh = fis.rce.build_hexahedral_mesh(nx, ny, nz, corners)
    raw_1f, raw_6f, vols = [], [], []
    scan_1f, scan_6f = [], []
    for beta in betas:
        A, M = fis.rce.assemble(mesh, beta)
        vals, vecs = fis.rce.solve_modes(A, M, modes)
        sets = fis.select_mode_sets(mesh, vals, vecs, beta, center_radius=center_radius)
        for label, idxs in sets.items():
            for mi in idxs:
                mode = vecs[:, mi]
                face_sigs = {f: fis.local_face_signature(mesh, mode, f, radius=face_radius, thickness=face_thickness) for f in FACES}
                face_grids = {f: fps.face_grid(mesh, mode, f) for f in FACES}
                for f in FACES:
                    raw_1f.append((fise.improved_face_vector(face_sigs[f]), label))
                    scan_1f.append((fps.scan_face_feature(face_grids[f], patch_sizes=(3,5)), label))
                raw_6f.append((fise.improved_sixface_vector(face_sigs), label))
                scan_6f.append((fps.scan_sixface_feature(face_grids, patch_sizes=(3,5)), label))
                vols.append((fise.improved_volume_vector(mesh, mode, center_radius=center_radius), label))
    return {
        'cube_raw_1F': raw_1f,
        'cube_raw_6F': raw_6f,
        'cube_scan_1F': scan_1f,
        'cube_scan_6F': scan_6f,
        'cube_volume': vols,
    }

def build_supercell_datasets(betas=(0.0,1.0,5.0), ncell=3, pts_per_cell=6, modes=20, radius=0.28, face_radius=0.24, face_thickness=0.08, topk=3):
    raw_1f, raw_6f, vols = [], [], []
    scan_1f, scan_6f = [], []
    for beta in betas:
        meshp, valsp, vecsp, rowsp, bestp, symp = pcs.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_prim = sorted(rowsp, key=lambda r: r['score'], reverse=True)[:topk]
        for r in top_prim:
            mi = r['mode'] - 1
            mode = vecsp[:, mi]
            fsigs = fpss.supercell_outer_face_sigs(meshp, mode, radius=face_radius, thickness=face_thickness)
            fgrids = {f: fps.face_grid(meshp, mode, f) for f in FACES}
            for f in FACES:
                raw_1f.append((fise.improved_face_vector(fsigs[f]), 'primitive_q'))
                scan_1f.append((fps.scan_face_feature(fgrids[f], patch_sizes=(3,5)), 'primitive_q'))
            raw_6f.append((fise.improved_sixface_vector(fsigs), 'primitive_q'))
            scan_6f.append((fps.scan_sixface_feature(fgrids, patch_sizes=(3,5)), 'primitive_q'))
            vols.append((fpss.primitive_volume_feature(r), 'primitive_q'))

        meshf, valsf, vecsf, rowsf, bestc, bestf, symf = fccr.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_fcc = sorted(rowsf, key=lambda r: r['fcc_score'], reverse=True)[:topk]
        for r in top_fcc:
            mi = r['mode'] - 1
            mode = vecsf[:, mi]
            fsigs = fpss.supercell_outer_face_sigs(meshf, mode, radius=face_radius, thickness=face_thickness)
            fgrids = {f: fps.face_grid(meshf, mode, f) for f in FACES}
            for f in FACES:
                raw_1f.append((fise.improved_face_vector(fsigs[f]), 'fcc_face'))
                scan_1f.append((fps.scan_face_feature(fgrids[f], patch_sizes=(3,5)), 'fcc_face'))
            raw_6f.append((fise.improved_sixface_vector(fsigs), 'fcc_face'))
            scan_6f.append((fps.scan_sixface_feature(fgrids, patch_sizes=(3,5)), 'fcc_face'))
            vols.append((fpss.fcc_volume_feature(r), 'fcc_face'))
    return {
        'supercell_raw_1F': raw_1f,
        'supercell_raw_6F': raw_6f,
        'supercell_scan_1F': scan_1f,
        'supercell_scan_6F': scan_6f,
        'supercell_volume': vols,
    }

def run_all(n_perm=300, seed=0):
    datasets = {}
    datasets.update(build_cube_datasets())
    datasets.update(build_supercell_datasets())
    rows = []
    for name, data in datasets.items():
        X = [x for x,_ in data]
        y = [lab for _,lab in data]
        res = permutation_test(X, y, n_perm=n_perm, seed=seed)
        res['dataset'] = name
        rows.append(res)
    return rows

def main():
    ap = argparse.ArgumentParser(description='Permutation/label-shuffle tests for current Helmholtz-Robin classification tasks.')
    ap.add_argument('--n-perm', type=int, default=300)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--csv', type=str, default='/mnt/data/label_shuffle_permutation_tests.csv')
    args = ap.parse_args()

    rows = run_all(n_perm=args.n_perm, seed=args.seed)
    out = Path(args.csv)
    with out.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['dataset','A_real','A_perm_mean','A_perm_std','A_perm_q95','A_perm_q99','p_ge','z_score','n_samples','n_classes'])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print('Label-shuffle permutation tests')
    for r in rows:
        print(f"{r['dataset']}: A_real={r['A_real']:.4f}, shuffle={r['A_perm_mean']:.4f}±{r['A_perm_std']:.4f}, q95={r['A_perm_q95']:.4f}, p_ge={r['p_ge']:.4f}, z={r['z_score']:.2f}")
    print(f'CSV written to {out}')

if __name__ == '__main__':
    main()
