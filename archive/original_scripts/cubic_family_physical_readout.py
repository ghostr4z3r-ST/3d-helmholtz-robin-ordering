import argparse, importlib.util, csv
from pathlib import Path
import pandas as pd


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

PCS = load_module('pcs', Path('/mnt/data/primitive_cubic_supercell_ordered.py'))
FCC = load_module('fcc', Path('/mnt/data/fcc_supercell_readout.py'))
BCC = load_module('bcc', Path('/mnt/data/bcc_supercell_readout.py'))


def pcf(center_xyz: float, face_xyz: float) -> float:
    return (center_xyz - face_xyz) / (center_xyz + face_xyz + 1e-15)


def pbc(body_xyz: float, corner_xyz: float) -> float:
    return (body_xyz - corner_xyz) / (body_xyz + corner_xyz + 1e-15)


def analyze_beta(beta: float, ncell: int = 3, pts_per_cell: int = 6, modes: int = 20, radius: float = 0.28):
    # Primitive cubic ordered readout
    _, _, _, prow_rows, pbest, _ = PCS.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)

    # FCC and BCC readouts
    _, _, _, frows, best_center, best_fcc, _ = FCC.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)
    _, _, _, brows, best_body, best_bcc, _ = BCC.analyze(beta, ncell, pts_per_cell, modes, radius, skip_ground=True)

    # Evaluate primitive mode on center-vs-face readout using the same mode index
    p_mode = pbest['mode']
    p_on_fcc = next(r for r in frows if r['mode'] == p_mode)

    row = {
        'beta': beta,
        'primitive_mode': p_mode,
        'primitive_lambda': pbest['lambda'],
        'primitive_center_xyz': p_on_fcc['center_xyz_mean_abs'],
        'primitive_face_xyz': p_on_fcc['fcc_face_mean'],
        'P_cf_primitive': pcf(p_on_fcc['center_xyz_mean_abs'], p_on_fcc['fcc_face_mean']),
        'primitive_q_label': pbest['best_q_label'],
        'primitive_q': str(pbest['best_q']),
        'primitive_q_amp': pbest['best_q_amp'],
        'primitive_uniformity': pbest['xyz_uniformity'],
        'primitive_same_sign': pbest['xyz_same_sign'],
        'primitive_score': pbest['score'],

        'fcc_mode': best_fcc['mode'],
        'fcc_lambda': best_fcc['lambda'],
        'fcc_center_xyz': best_fcc['center_xyz_mean_abs'],
        'fcc_face_xyz': best_fcc['fcc_face_mean'],
        'P_cf_fcc': pcf(best_fcc['center_xyz_mean_abs'], best_fcc['fcc_face_mean']),
        'fcc_q_xy': best_fcc['xy_best_q_label'],
        'fcc_q_xz': best_fcc['xz_best_q_label'],
        'fcc_q_yz': best_fcc['yz_best_q_label'],
        'fcc_uniformity': best_fcc['fcc_face_uniformity'],
        'fcc_score': best_fcc['fcc_score'],

        'bcc_mode': best_bcc['mode'],
        'bcc_lambda': best_bcc['lambda'],
        'bcc_body_xyz': best_bcc['body_xyz_mean_abs'],
        'bcc_corner_xyz': best_bcc['corner_xyz_mean_abs'],
        'P_bc_bcc': pbc(best_bcc['body_xyz_mean_abs'], best_bcc['corner_xyz_mean_abs']),
        'bcc_q_body': best_bcc['body_best_q_label'],
        'bcc_q_corner': best_bcc['corner_best_q_label'],
        'bcc_body_uniformity': best_bcc['body_xyz_uniformity'],
        'bcc_corner_uniformity': best_bcc['corner_xyz_uniformity'],
        'bcc_score': best_bcc['bcc_score'],
    }
    return row


def main():
    ap = argparse.ArgumentParser(description='First physical proxy readout for primitive cubic vs bcc vs fcc.')
    ap.add_argument('--betas', type=str, default='0,1,3,4.5,5,10')
    ap.add_argument('--ncell', type=int, default=3)
    ap.add_argument('--pts-per-cell', type=int, default=6)
    ap.add_argument('--modes', type=int, default=20)
    ap.add_argument('--radius', type=float, default=0.28)
    ap.add_argument('--csv', type=str, default='')
    args = ap.parse_args()

    betas = [float(x) for x in args.betas.split(',') if x.strip()]
    rows = []
    for beta in betas:
        rows.append(analyze_beta(beta, args.ncell, args.pts_per_cell, args.modes, args.radius))

    df = pd.DataFrame(rows)
    with pd.option_context('display.max_columns', None, 'display.width', 220):
        print(df)

    if args.csv:
        df.to_csv(args.csv, index=False)


if __name__ == '__main__':
    main()
