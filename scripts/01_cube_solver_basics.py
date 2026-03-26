from __future__ import annotations

import argparse
import csv
from pathlib import Path

from paper1.cube_fem import robin_3d_matrices, solve_robin_cube, symmetry_error, convergence_table, beta_scan


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description='Core cube solver demo for the 3D Helmholtz-Robin eigenproblem.')
    ap.add_argument('--n', type=int, default=10)
    ap.add_argument('--L', type=float, default=1.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--modes', type=int, default=8)
    ap.add_argument('--tol', type=float, default=1e-9)
    ap.add_argument('--convergence', nargs='*', type=int, default=[6, 8, 10, 12])
    ap.add_argument('--beta-scan', nargs='*', type=float, dest='beta_scan_vals', default=[0.0, 0.1, 1.0, 10.0])
    ap.add_argument('--csv-dir', type=str, default='')
    args = ap.parse_args()

    A, _ = robin_3d_matrices(n=args.n, L=args.L, beta=args.beta)
    res = solve_robin_cube(n=args.n, L=args.L, beta=args.beta, modes=args.modes, tol=args.tol)
    print(f'Matrix symmetry ||A-A^T||_F: {symmetry_error(A):.3e}')
    print('Lowest eigenvalues:')
    for i, (lam, kval) in enumerate(zip(res.eigenvalues, res.k_values), start=1):
        print(f'  mode {i:2d}: lambda={lam:.12f}, k={kval:.12f}')

    conv = convergence_table(args.convergence, L=args.L, beta=args.beta, modes=min(args.modes, 6), tol=args.tol)
    print('\nConvergence table:')
    for n, vals in conv:
        print(f'  n={n}: ' + ', '.join(f'{v:.8f}' for v in vals[:6]))

    scan = beta_scan(args.beta_scan_vals, n=args.n, L=args.L, modes=min(args.modes, 6), tol=args.tol)
    print('\nBeta scan:')
    for beta, vals in scan:
        print(f'  beta={beta}: ' + ', '.join(f'{v:.8f}' for v in vals[:6]))

    if args.csv_dir:
        outdir = Path(args.csv_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        eig_rows = [{'mode': i + 1, 'lambda': float(lam), 'k': float(k)} for i, (lam, k) in enumerate(zip(res.eigenvalues, res.k_values))]
        write_rows(outdir / 'cube_eigenvalues.csv', eig_rows)
        conv_rows = []
        for n, vals in conv:
            row = {'n': n}
            for i, v in enumerate(vals[:6], start=1):
                row[f'lambda_{i}'] = float(v)
            conv_rows.append(row)
        write_rows(outdir / 'cube_convergence.csv', conv_rows)
        scan_rows = []
        for beta, vals in scan:
            row = {'beta': beta}
            for i, v in enumerate(vals[:6], start=1):
                row[f'lambda_{i}'] = float(v)
            scan_rows.append(row)
        write_rows(outdir / 'cube_beta_scan.csv', scan_rows)


if __name__ == '__main__':
    main()
