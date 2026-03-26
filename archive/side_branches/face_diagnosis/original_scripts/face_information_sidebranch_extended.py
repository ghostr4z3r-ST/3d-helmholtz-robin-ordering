import argparse
import csv
import importlib.util
from pathlib import Path
from collections import defaultdict
import numpy as np

# imports
spec = importlib.util.spec_from_file_location('fis', '/mnt/data/face_information_sidebranch.py')
fis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fis)

spec2 = importlib.util.spec_from_file_location('pcs', '/mnt/data/primitive_cubic_supercell_ordered.py')
pcs = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(pcs)

spec3 = importlib.util.spec_from_file_location('fccr', '/mnt/data/fcc_supercell_readout.py')
fccr = importlib.util.module_from_spec(spec3)
spec3.loader.exec_module(fccr)

FACES = ['x0','x1','y0','y1','z0','z1']


def improved_face_vector(sig):
    coeffs = sig['coeffs']
    vals = np.asarray(sig['values4'], dtype=float)
    valsc = vals - vals.mean()
    energy = float(np.linalg.norm(vals))
    absmean = float(np.mean(np.abs(vals)))
    sign_imb = float(abs(np.sum(np.sign(vals + 1e-15))) / 4.0)
    spread = float(np.std(np.abs(vals)))
    # raw and abs coefficients both retained
    return np.array([
        energy, absmean, sign_imb, spread,
        coeffs['const'], coeffs['a'], coeffs['b'], coeffs['ab'],
        abs(coeffs['const']), abs(coeffs['a']), abs(coeffs['b']), abs(coeffs['ab'])
    ], dtype=float)


def improved_sixface_vector(face_sigs):
    vecs = [improved_face_vector(face_sigs[f]) for f in FACES]
    arr = np.vstack(vecs)
    mean = arr.mean(axis=0)
    std = arr.std(axis=0)
    # opposite-face correlations on coeff vectors (excluding scalar summaries)
    pairs = [('x0','x1'), ('y0','y1'), ('z0','z1')]
    cors = []
    for a,b in pairs:
        va = improved_face_vector(face_sigs[a])[4:]
        vb = improved_face_vector(face_sigs[b])[4:]
        na = np.linalg.norm(va); nb = np.linalg.norm(vb)
        cors.append(float(va @ vb / (na * nb + 1e-15)))
    return np.concatenate([mean, std, np.array(cors, dtype=float)])


def improved_volume_vector(mesh, vec, center_radius=0.22):
    coeffs = fis.center_signature(mesh, vec, radius=center_radius)
    raw = np.array([coeffs[k] for k in ['x','y','z','xy','xz','yz','xyz']], dtype=float)
    return np.concatenate([raw, np.abs(raw), [np.linalg.norm(raw)]])


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
        correct += (pred == s['label'])
    return correct / len(samples)


def cube_experiment(betas=(0.0,1.0,3.0,5.0), nx=17, ny=17, nz=17, modes=14, face_radius=0.24, face_thickness=0.08, center_radius=0.22):
    corners = fis.rce.make_corners('cube')
    mesh = fis.rce.build_hexahedral_mesh(nx, ny, nz, corners)
    samples_1f, samples_6f, samples_v = [], [], []
    rows = []
    for beta in betas:
        A, M = fis.rce.assemble(mesh, beta)
        vals, vecs = fis.rce.solve_modes(A, M, modes)
        sets = fis.select_mode_sets(mesh, vals, vecs, beta, center_radius=center_radius)
        # take all selected modes (more modes per class than first test)
        for label, idxs in sets.items():
            for mi in idxs:
                mode = vecs[:, mi]
                face_sigs = {f: fis.local_face_signature(mesh, mode, f, radius=face_radius, thickness=face_thickness) for f in FACES}
                for f in FACES:
                    samples_1f.append({'label': label, 'feature': improved_face_vector(face_sigs[f])})
                samples_6f.append({'label': label, 'feature': improved_sixface_vector(face_sigs)})
                samples_v.append({'label': label, 'feature': improved_volume_vector(mesh, mode, center_radius=center_radius)})
                rows.append({
                    'experiment': 'cube', 'beta': beta, 'label': label, 'mode': mi+1, 'lambda': float(vals[mi]),
                    'center_top': fis.top_abs_coeff_name(fis.center_signature(mesh, mode, radius=center_radius)),
                    'face_pattern': '|'.join(fis.top_abs_coeff_name(face_sigs[f]['coeffs']) for f in FACES)
                })
    A1 = nearest_centroid_cv(samples_1f)
    A6 = nearest_centroid_cv(samples_6f)
    AV = nearest_centroid_cv(samples_v)
    return rows, {'A_1F': A1, 'A_6F': A6, 'A_V': AV, 'D_1F': 1 - A1/AV if AV>0 else np.nan, 'D_6F': 1 - A6/AV if AV>0 else np.nan,
                  'n_face_samples': len(samples_1f), 'n_multiface_samples': len(samples_6f), 'n_volume_samples': len(samples_v)}


def supercell_outer_face_sigs(mesh, modevec, radius=0.24, thickness=0.08):
    return {f: fis.local_face_signature(mesh, modevec, f, radius=radius, thickness=thickness) for f in FACES}


def primitive_volume_feature(row):
    return np.array([
        row['xyz_mean_abs'], row['xyz_uniformity'], row['xyz_same_sign'], row['best_q_amp'],
        row['corr_x'], row['corr_y'], row['corr_z'], row['axis_mean_abs'], row['pair_mean_abs'],
        row['dominance_over_axis'], row['dominance_over_pair']
    ], dtype=float)


def fcc_volume_feature(row):
    # capture center-vs-face and overall face organization
    pcf = (row['center_xyz_mean_abs'] - row['fcc_face_mean']) / (row['center_xyz_mean_abs'] + row['fcc_face_mean'] + 1e-15)
    return np.array([
        row['center_xyz_mean_abs'], row['fcc_face_mean'], pcf,
        row['center_best_q_amp'], row['fcc_face_qamp'],
        row['xy_xyz_mean_abs'], row['xz_xyz_mean_abs'], row['yz_xyz_mean_abs'],
        row['fcc_face_uniformity'], row['center_xyz_uniformity'],
        (row['xy_best_q_amp'] + row['xz_best_q_amp'] + row['yz_best_q_amp']) / 3.0
    ], dtype=float)


def supercell_experiment(betas=(0.0,1.0,5.0), ncell=3, pts_per_cell=6, modes=20, radius=0.28):
    samples_1f, samples_6f, samples_v = [], [], []
    rows = []
    for beta in betas:
        # primitive cubic top 3 by score
        meshp, valsp, vecsp, rowsp, bestp, symp = pcs.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_prim = sorted(rowsp, key=lambda r: r['score'], reverse=True)[:3]
        for r in top_prim:
            mi = r['mode'] - 1
            fsigs = supercell_outer_face_sigs(meshp, vecsp[:, mi], radius=0.24, thickness=0.08)
            for f in FACES:
                samples_1f.append({'label': 'primitive_q', 'feature': improved_face_vector(fsigs[f])})
            samples_6f.append({'label': 'primitive_q', 'feature': improved_sixface_vector(fsigs)})
            samples_v.append({'label': 'primitive_q', 'feature': primitive_volume_feature(r)})
            rows.append({'experiment':'supercell','family':'primitive_q','beta':beta,'mode':r['mode'],'lambda':r['lambda'],'best_q':r['best_q_label']})
        # fcc top 3 by fcc_score
        meshf, valsf, vecsf, rowsf, bestc, bestf, symf = fccr.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_fcc = sorted(rowsf, key=lambda r: r['fcc_score'], reverse=True)[:3]
        for r in top_fcc:
            mi = r['mode'] - 1
            fsigs = supercell_outer_face_sigs(meshf, vecsf[:, mi], radius=0.24, thickness=0.08)
            for f in FACES:
                samples_1f.append({'label': 'fcc_face', 'feature': improved_face_vector(fsigs[f])})
            samples_6f.append({'label': 'fcc_face', 'feature': improved_sixface_vector(fsigs)})
            samples_v.append({'label': 'fcc_face', 'feature': fcc_volume_feature(r)})
            rows.append({'experiment':'supercell','family':'fcc_face','beta':beta,'mode':r['mode'],'lambda':r['lambda'],'best_q':r['xy_best_q_label']})
    A1 = nearest_centroid_cv(samples_1f)
    A6 = nearest_centroid_cv(samples_6f)
    AV = nearest_centroid_cv(samples_v)
    return rows, {'A_1F': A1, 'A_6F': A6, 'A_V': AV, 'D_1F': 1 - A1/AV if AV>0 else np.nan, 'D_6F': 1 - A6/AV if AV>0 else np.nan,
                  'n_face_samples': len(samples_1f), 'n_multiface_samples': len(samples_6f), 'n_volume_samples': len(samples_v)}


def main():
    ap = argparse.ArgumentParser(description='Extended sidebranch test for face information deficit using improved features and supercells.')
    ap.add_argument('--csv', type=str, default='/mnt/data/face_information_extended_test.csv')
    args = ap.parse_args()

    cube_rows, cube_metrics = cube_experiment()
    super_rows, super_metrics = supercell_experiment()

    out = Path(args.csv)
    with out.open('w', newline='') as f:
        fieldnames = ['section','experiment','beta','label_or_family','mode','lambda','note','value']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in cube_rows:
            w.writerow({'section':'cube_rows','experiment':r['experiment'],'beta':r['beta'],'label_or_family':r['label'],'mode':r['mode'],'lambda':r['lambda'],'note':f"center={r['center_top']}; faces={r['face_pattern']}",'value':''})
        for k,v in cube_metrics.items():
            w.writerow({'section':'cube_metrics','experiment':'cube','beta':'','label_or_family':k,'mode':'','lambda':'','note':'','value':v})
        for r in super_rows:
            w.writerow({'section':'super_rows','experiment':r['experiment'],'beta':r['beta'],'label_or_family':r['family'],'mode':r['mode'],'lambda':r['lambda'],'note':f"best_q={r['best_q']}",'value':''})
        for k,v in super_metrics.items():
            w.writerow({'section':'super_metrics','experiment':'supercell','beta':'','label_or_family':k,'mode':'','lambda':'','note':'','value':v})

    print('Extended face-information sidebranch test')
    print('Cube metrics:')
    for k,v in cube_metrics.items():
        print(f'  {k} = {v}')
    print('Supercell primitive-vs-fcc metrics:')
    for k,v in super_metrics.items():
        print(f'  {k} = {v}')
    print(f'CSV written to {out}')

if __name__ == '__main__':
    main()
