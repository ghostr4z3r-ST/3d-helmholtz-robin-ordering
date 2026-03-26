import argparse, csv, importlib.util
from pathlib import Path
import numpy as np

# import prior modules
spec0 = importlib.util.spec_from_file_location('fps', '/mnt/data/face_patch_scan_sidebranch.py')
fps = importlib.util.module_from_spec(spec0); spec0.loader.exec_module(fps)

spec1 = importlib.util.spec_from_file_location('fise', '/mnt/data/face_information_sidebranch_extended.py')
fise = importlib.util.module_from_spec(spec1); spec1.loader.exec_module(fise)

spec2 = importlib.util.spec_from_file_location('pcs', '/mnt/data/primitive_cubic_supercell_ordered.py')
pcs = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(pcs)

spec3 = importlib.util.spec_from_file_location('fccr', '/mnt/data/fcc_supercell_readout.py')
fccr = importlib.util.module_from_spec(spec3); spec3.loader.exec_module(fccr)

FACES = ['x0','x1','y0','y1','z0','z1']


def supercell_outer_face_sigs(mesh, modevec, radius=0.24, thickness=0.08):
    return {f: fise.fis.local_face_signature(mesh, modevec, f, radius=radius, thickness=thickness) for f in FACES}


def primitive_volume_feature(row):
    return np.array([
        row['xyz_mean_abs'], row['xyz_uniformity'], row['xyz_same_sign'], row['best_q_amp'],
        row['corr_x'], row['corr_y'], row['corr_z'], row['axis_mean_abs'], row['pair_mean_abs'],
        row['dominance_over_axis'], row['dominance_over_pair']
    ], dtype=float)


def fcc_volume_feature(row):
    pcf = (row['center_xyz_mean_abs'] - row['fcc_face_mean']) / (row['center_xyz_mean_abs'] + row['fcc_face_mean'] + 1e-15)
    return np.array([
        row['center_xyz_mean_abs'], row['fcc_face_mean'], pcf,
        row['center_best_q_amp'], row['fcc_face_qamp'],
        row['xy_xyz_mean_abs'], row['xz_xyz_mean_abs'], row['yz_xyz_mean_abs'],
        row['fcc_face_uniformity'], row['center_xyz_uniformity'],
        (row['xy_best_q_amp'] + row['xz_best_q_amp'] + row['yz_best_q_amp']) / 3.0
    ], dtype=float)


def run_supercell_scan(betas=(0.0,1.0,5.0), ncell=3, pts_per_cell=6, modes=20, radius=0.28, face_radius=0.24, face_thickness=0.08, patch_sizes=(3,5), topk=3):
    raw_1f, raw_6f, scan_1f, scan_6f, vols = [], [], [], [], []
    rows = []
    for beta in betas:
        meshp, valsp, vecsp, rowsp, bestp, symp = pcs.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_prim = sorted(rowsp, key=lambda r: r['score'], reverse=True)[:topk]
        for r in top_prim:
            mi = r['mode'] - 1
            mode = vecsp[:, mi]
            fsigs = supercell_outer_face_sigs(meshp, mode, radius=face_radius, thickness=face_thickness)
            fgrids = {f: fps.face_grid(meshp, mode, f) for f in FACES}
            for f in FACES:
                raw_1f.append({'label': 'primitive_q', 'feature': fise.improved_face_vector(fsigs[f])})
                scan_1f.append({'label': 'primitive_q', 'feature': fps.scan_face_feature(fgrids[f], patch_sizes=patch_sizes)})
            raw_6f.append({'label': 'primitive_q', 'feature': fise.improved_sixface_vector(fsigs)})
            scan_6f.append({'label': 'primitive_q', 'feature': fps.scan_sixface_feature(fgrids, patch_sizes=patch_sizes)})
            vols.append({'label': 'primitive_q', 'feature': primitive_volume_feature(r)})
            rows.append({'family':'primitive_q','beta':beta,'mode':r['mode'],'lambda':r['lambda'],'note':f"best_q={r['best_q_label']}"})

        meshf, valsf, vecsf, rowsf, bestc, bestf, symf = fccr.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
        top_fcc = sorted(rowsf, key=lambda r: r['fcc_score'], reverse=True)[:topk]
        for r in top_fcc:
            mi = r['mode'] - 1
            mode = vecsf[:, mi]
            fsigs = supercell_outer_face_sigs(meshf, mode, radius=face_radius, thickness=face_thickness)
            fgrids = {f: fps.face_grid(meshf, mode, f) for f in FACES}
            for f in FACES:
                raw_1f.append({'label': 'fcc_face', 'feature': fise.improved_face_vector(fsigs[f])})
                scan_1f.append({'label': 'fcc_face', 'feature': fps.scan_face_feature(fgrids[f], patch_sizes=patch_sizes)})
            raw_6f.append({'label': 'fcc_face', 'feature': fise.improved_sixface_vector(fsigs)})
            scan_6f.append({'label': 'fcc_face', 'feature': fps.scan_sixface_feature(fgrids, patch_sizes=patch_sizes)})
            vols.append({'label': 'fcc_face', 'feature': fcc_volume_feature(r)})
            rows.append({'family':'fcc_face','beta':beta,'mode':r['mode'],'lambda':r['lambda'],'note':f"xy_q={r['xy_best_q_label']}"})

    metrics = {}
    metrics['A_1F'] = fise.nearest_centroid_cv(raw_1f)
    metrics['A_6F'] = fise.nearest_centroid_cv(raw_6f)
    metrics['A_1F_scan'] = fise.nearest_centroid_cv(scan_1f)
    metrics['A_6F_scan'] = fise.nearest_centroid_cv(scan_6f)
    metrics['A_V'] = fise.nearest_centroid_cv(vols)
    for k in ['A_1F', 'A_6F', 'A_1F_scan', 'A_6F_scan']:
        metrics['D_' + k[2:]] = 1.0 - metrics[k] / metrics['A_V'] if metrics['A_V'] > 0 else float('nan')
    metrics['G_1F'] = metrics['A_1F_scan'] - metrics['A_1F']
    metrics['G_6F'] = metrics['A_6F_scan'] - metrics['A_6F']
    metrics['n_face_samples'] = len(raw_1f)
    metrics['n_multiface_samples'] = len(raw_6f)
    metrics['n_volume_samples'] = len(vols)
    return rows, metrics


def main():
    ap = argparse.ArgumentParser(description='Patch-scan sidebranch on primitive cubic supercell vs fcc.')
    ap.add_argument('--betas', type=str, default='0,1,5')
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--patch-sizes', type=str, default='3,5')
    ap.add_argument('--topk', type=int, default=3)
    ap.add_argument('--csv', type=str, default='/mnt/data/face_patch_scan_supercell_test.csv')
    args = ap.parse_args()

    betas = tuple(float(x) for x in args.betas.split(','))
    patch_sizes = tuple(int(x) for x in args.patch_sizes.split(','))
    rows, metrics = run_supercell_scan(betas=betas, ncell=args.ncell, pts_per_cell=args.pts_per_cell, modes=args.modes, patch_sizes=patch_sizes, topk=args.topk)

    out = Path(args.csv)
    with out.open('w', newline='') as f:
        fieldnames = ['section','family','beta','mode','lambda','note','metric','value']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({'section':'rows', **r, 'metric':'', 'value':''})
        for k,v in metrics.items():
            w.writerow({'section':'metrics','family':'','beta':'','mode':'','lambda':'','note':'','metric':k,'value':v})

    print('Supercell patch-scan sidebranch')
    for k,v in metrics.items():
        print(f'  {k} = {v}')
    print(f'CSV written to {out}')

if __name__ == '__main__':
    main()
