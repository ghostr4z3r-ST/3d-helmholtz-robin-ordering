
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
pcs = load('pcs', '/mnt/data/primitive_cubic_supercell_ordered.py')
fccr = load('fccr', '/mnt/data/fcc_supercell_readout.py')

FACES = ['x0','x1','y0','y1','z0','z1']

def nc_acc(X, y):
    X = [np.asarray(v, float) for v in X]
    y = list(y)
    n = len(y)
    c = 0
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        labs = sorted(set(y[j] for j in tr))
        cents = {}
        for lab in labs:
            cents[lab] = np.vstack([X[j] for j in tr if y[j]==lab]).mean(axis=0)
        pred = min(labs, key=lambda lab: np.linalg.norm(X[i]-cents[lab]))
        c += int(pred == y[i])
    return c / n

def groups_for(name, dim):
    if name.endswith('raw_1F'):
        return {'summaries': range(0,4), 'signed_coeffs': range(4,8), 'abs_coeffs': range(8,12)}
    if name.endswith('raw_6F'):
        return {'face_means': range(0,12), 'face_stds': range(12,24), 'opp_corr': range(24,27)}
    if name == 'cube_volume':
        return {'signed_axis': range(0,3), 'signed_pair3d': range(3,7), 'abs_axis': range(7,10), 'abs_pair3d': range(10,14), 'norm':[14]}
    if name.endswith('scan_1F'):
        return {'patch_mean': range(0,17), 'patch_std': range(17,34), 'patch_q25': range(34,51), 'patch_q75': range(51,68)}
    if name == 'supercell_volume':
        return {'core_amp': range(0,3), 'q_or_balance': range(3,5), 'subcomponents': range(5,8), 'uniformity_tail': range(8,11)}
    raise ValueError((name, dim))

def eval_ablation(name, X, y):
    X = [np.asarray(v,float) for v in X]
    full = nc_acc(X, y)
    d = len(X[0])
    rows = [{'dataset': name, 'kind': 'full', 'group':'ALL', 'n_features':d, 'accuracy': full, 'delta_from_full': 0.0}]
    groups = groups_for(name, d)
    for g, idx in groups.items():
        idx = list(idx)
        only = [x[idx] for x in X]
        a_only = nc_acc(only, y)
        keep = [i for i in range(d) if i not in set(idx)]
        a_drop = nc_acc([x[keep] for x in X], y) if keep else float('nan')
        rows.append({'dataset': name, 'kind': 'only', 'group': g, 'n_features': len(idx), 'accuracy': a_only, 'delta_from_full': a_only-full})
        rows.append({'dataset': name, 'kind': 'drop', 'group': g, 'n_features': len(keep), 'accuracy': a_drop, 'delta_from_full': a_drop-full})
    return rows

def build_cube_small():
    corners = fis.rce.make_corners('cube')
    mesh = fis.rce.build_hexahedral_mesh(13,13,13,corners)
    X = {'cube_raw_1F': [], 'cube_raw_6F': [], 'cube_scan_1F': [], 'cube_volume': []}
    y = {k: [] for k in X}
    for beta in (0.0,1.0,5.0):
        A, M = fis.rce.assemble(mesh, beta)
        vals, vecs = fis.rce.solve_modes(A, M, 12)
        sets = fis.select_mode_sets(mesh, vals, vecs, beta, center_radius=0.22)
        for lab, idxs in sets.items():
            for mi in idxs:
                mode = vecs[:, mi]
                face_sigs = {f: fis.local_face_signature(mesh, mode, f, radius=0.24, thickness=0.08) for f in FACES}
                face_grids = {f: fps.face_grid(mesh, mode, f) for f in FACES}
                for f in FACES:
                    X['cube_raw_1F'].append(fise.improved_face_vector(face_sigs[f])); y['cube_raw_1F'].append(lab)
                    X['cube_scan_1F'].append(fps.scan_face_feature(face_grids[f], patch_sizes=(3,5))); y['cube_scan_1F'].append(lab)
                X['cube_raw_6F'].append(fise.improved_sixface_vector(face_sigs)); y['cube_raw_6F'].append(lab)
                X['cube_volume'].append(fise.improved_volume_vector(mesh, mode, center_radius=0.22)); y['cube_volume'].append(lab)
    return X, y

def build_supercell_small():
    X = {'supercell_raw_1F': [], 'supercell_scan_1F': [], 'supercell_volume': []}
    y = {k: [] for k in X}
    for beta in (0.0,1.0,5.0):
        meshp, valsp, vecsp, rowsp, bestp, symp = pcs.analyze(beta, 3, 5, 16, 0.28, skip_ground=True)
        top_prim = sorted(rowsp, key=lambda r: r['score'], reverse=True)[:2]
        for r in top_prim:
            mode = vecsp[:, r['mode']-1]
            fsigs = {f: fis.local_face_signature(meshp, mode, f, radius=0.24, thickness=0.08) for f in FACES}
            fgrids = {f: fps.face_grid(meshp, mode, f) for f in FACES}
            for f in FACES:
                X['supercell_raw_1F'].append(fise.improved_face_vector(fsigs[f])); y['supercell_raw_1F'].append('primitive_q')
                X['supercell_scan_1F'].append(fps.scan_face_feature(fgrids[f], patch_sizes=(3,5))); y['supercell_scan_1F'].append('primitive_q')
            X['supercell_volume'].append(np.array([
                r['xyz_mean_abs'], r['xyz_uniformity'], r['xyz_same_sign'], r['best_q_amp'],
                r['corr_x'], r['corr_y'], r['corr_z'], r['axis_mean_abs'], r['pair_mean_abs'],
                r['dominance_over_axis'], r['dominance_over_pair']], float)); y['supercell_volume'].append('primitive_q')
        meshf, valsf, vecsf, rowsf, bestc, bestf, symf = fccr.analyze(beta, 3, 5, 16, 0.28, skip_ground=True)
        top_fcc = sorted(rowsf, key=lambda r: r['fcc_score'], reverse=True)[:2]
        for r in top_fcc:
            mode = vecsf[:, r['mode']-1]
            fsigs = {f: fis.local_face_signature(meshf, mode, f, radius=0.24, thickness=0.08) for f in FACES}
            fgrids = {f: fps.face_grid(meshf, mode, f) for f in FACES}
            for f in FACES:
                X['supercell_raw_1F'].append(fise.improved_face_vector(fsigs[f])); y['supercell_raw_1F'].append('fcc_face')
                X['supercell_scan_1F'].append(fps.scan_face_feature(fgrids[f], patch_sizes=(3,5))); y['supercell_scan_1F'].append('fcc_face')
            pcf = (r['center_xyz_mean_abs'] - r['fcc_face_mean']) / (r['center_xyz_mean_abs'] + r['fcc_face_mean'] + 1e-15)
            X['supercell_volume'].append(np.array([
                r['center_xyz_mean_abs'], r['fcc_face_mean'], pcf,
                r['center_best_q_amp'], r['fcc_face_qamp'],
                r['xy_xyz_mean_abs'], r['xz_xyz_mean_abs'], r['yz_xyz_mean_abs'],
                r['fcc_face_uniformity'], r['center_xyz_uniformity'],
                (r['xy_best_q_amp'] + r['xz_best_q_amp'] + r['yz_best_q_amp'])/3.0], float)); y['supercell_volume'].append('fcc_face')
    return X, y

def summarize(detail_rows):
    out = []
    by = {}
    for r in detail_rows:
        by.setdefault(r['dataset'], []).append(r)
    for ds, rs in by.items():
        full = next(r['accuracy'] for r in rs if r['kind']=='full')
        drops = [r for r in rs if r['kind']=='drop' and np.isfinite(r['accuracy'])]
        onlys = [r for r in rs if r['kind']=='only']
        worst_drop = min(drops, key=lambda r: r['delta_from_full'])
        best_only = max(onlys, key=lambda r: r['accuracy'])
        out.append({
            'dataset': ds,
            'full_accuracy': full,
            'most_informative_group': worst_drop['group'],
            'accuracy_without_group': worst_drop['accuracy'],
            'drop_delta': worst_drop['delta_from_full'],
            'best_single_group': best_only['group'],
            'best_single_accuracy': best_only['accuracy'],
            'single_delta': best_only['delta_from_full'],
        })
    return out

def main():
    detail=[]
    Xc, yc = build_cube_small()
    for name in Xc:
        detail.extend(eval_ablation(name, Xc[name], yc[name]))
    Xs, ys = build_supercell_small()
    for name in Xs:
        detail.extend(eval_ablation(name, Xs[name], ys[name]))
    summary = summarize(detail)
    with open('/mnt/data/feature_ablation_tests_reduced.csv','w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=['dataset','kind','group','n_features','accuracy','delta_from_full'])
        w.writeheader(); w.writerows(detail)
    with open('/mnt/data/feature_ablation_summary_reduced.csv','w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=['dataset','full_accuracy','most_informative_group','accuracy_without_group','drop_delta','best_single_group','best_single_accuracy','single_delta'])
        w.writeheader(); w.writerows(summary)
    print('done')
    for r in summary:
        print(r)

if __name__=='__main__':
    main()
