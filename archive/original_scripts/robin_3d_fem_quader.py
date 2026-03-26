from __future__ import annotations

"""
3D Robin eigenvalue problem on a rectangular box (quader) using tensor-product
linear finite elements.

We solve
    -Δu = λ u  in Ω = [0, Lx] x [0, Ly] x [0, Lz]
with Robin boundary condition
    ∂_n u + beta * u = 0  on ∂Ω.

Weak form:
    ∫_Ω ∇u·∇v dΩ + beta ∫_{∂Ω} u v dS = λ ∫_Ω u v dΩ

This yields the symmetric generalized eigenproblem
    A u = λ M u
with tensor-product FEM matrices.
"""

import argparse
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class EigenResult:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    k_values: np.ndarray


def fem1d_matrices(n: int, L: float) -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix]:
    if n < 2:
        raise ValueError("n must be at least 2.")
    h = L / (n - 1)

    main_k = np.full(n, 2.0 / h)
    off_k = np.full(n - 1, -1.0 / h)
    K = sp.diags([off_k, main_k, off_k], offsets=[-1, 0, 1], format="lil")
    K[0, 0] = 1.0 / h
    K[-1, -1] = 1.0 / h
    K = K.tocsr()

    main_m = np.full(n, 4.0 * h / 6.0)
    off_m = np.full(n - 1, 1.0 * h / 6.0)
    M = sp.diags([off_m, main_m, off_m], offsets=[-1, 0, 1], format="lil")
    M[0, 0] = h / 3.0
    M[-1, -1] = h / 3.0
    M = M.tocsr()

    bdiag = np.zeros(n, dtype=float)
    bdiag[0] = 1.0
    bdiag[-1] = 1.0
    B = sp.diags(bdiag, 0, format="csr")

    return K, M, B


def robin_3d_quader_matrices(
    nx: int,
    ny: int,
    nz: int,
    Lx: float,
    Ly: float,
    Lz: float,
    beta: float,
) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    Kx, Mx, Bx = fem1d_matrices(nx, Lx)
    Ky, My, By = fem1d_matrices(ny, Ly)
    Kz, Mz, Bz = fem1d_matrices(nz, Lz)

    # Volume stiffness
    A_vol = (
        sp.kron(sp.kron(Kx, My, format="csr"), Mz, format="csr")
        + sp.kron(sp.kron(Mx, Ky, format="csr"), Mz, format="csr")
        + sp.kron(sp.kron(Mx, My, format="csr"), Kz, format="csr")
    )

    # Robin boundary terms on all six faces
    # x-faces: integrate over y-z surface -> My ⊗ Mz
    # y-faces: integrate over x-z surface -> Mx ⊗ Mz
    # z-faces: integrate over x-y surface -> Mx ⊗ My
    R = beta * (
        sp.kron(sp.kron(Bx, My, format="csr"), Mz, format="csr")
        + sp.kron(sp.kron(Mx, By, format="csr"), Mz, format="csr")
        + sp.kron(sp.kron(Mx, My, format="csr"), Bz, format="csr")
    )

    M = sp.kron(sp.kron(Mx, My, format="csr"), Mz, format="csr")
    return (A_vol + R).tocsr(), M.tocsr()


def solve_robin_quader(
    nx: int,
    ny: int,
    nz: int,
    Lx: float,
    Ly: float,
    Lz: float,
    beta: float,
    modes: int,
    tol: float = 1e-10,
) -> EigenResult:
    A, M = robin_3d_quader_matrices(nx, ny, nz, Lx, Ly, Lz, beta)
    ndof = A.shape[0]
    if modes >= ndof - 1:
        raise ValueError("modes must be smaller than ndof-1")
    vals, vecs = spla.eigsh(A, k=modes, M=M, sigma=0.0, which="LM", tol=tol)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    kvals = np.sqrt(np.clip(vals, 0.0, None))
    return EigenResult(vals, vecs, kvals)


def symmetry_error(A: sp.csr_matrix) -> float:
    D = A - A.T
    if D.nnz == 0:
        return 0.0
    return float(np.sqrt(np.sum(D.data * D.data)))


def convergence_table(
    ns: Iterable[int],
    Lx: float,
    Ly: float,
    Lz: float,
    beta: float,
    modes: int,
    tol: float,
) -> list[tuple[int, np.ndarray]]:
    out: list[tuple[int, np.ndarray]] = []
    for n in ns:
        res = solve_robin_quader(n, n, n, Lx, Ly, Lz, beta, modes, tol)
        out.append((n, res.eigenvalues.copy()))
    return out


def beta_scan(
    betas: Iterable[float],
    nx: int,
    ny: int,
    nz: int,
    Lx: float,
    Ly: float,
    Lz: float,
    modes: int,
    tol: float,
) -> list[tuple[float, np.ndarray]]:
    out: list[tuple[float, np.ndarray]] = []
    for beta in betas:
        res = solve_robin_quader(nx, ny, nz, Lx, Ly, Lz, beta, modes, tol)
        out.append((beta, res.eigenvalues.copy()))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Solve the 3D Robin eigenvalue problem on a rectangular box using tensor-product FEM.")
    parser.add_argument("--nx", type=int, default=10)
    parser.add_argument("--ny", type=int, default=10)
    parser.add_argument("--nz", type=int, default=10)
    parser.add_argument("--Lx", type=float, default=1.0)
    parser.add_argument("--Ly", type=float, default=1.2)
    parser.add_argument("--Lz", type=float, default=0.8)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--modes", type=int, default=12)
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--convergence", nargs="*", type=int)
    parser.add_argument("--beta-scan", nargs="*", type=float, dest="beta_scan_vals")
    args = parser.parse_args()

    A, _ = robin_3d_quader_matrices(args.nx, args.ny, args.nz, args.Lx, args.Ly, args.Lz, args.beta)
    print(f"Geometry: Lx={args.Lx}, Ly={args.Ly}, Lz={args.Lz}")
    print(f"Grid: nx={args.nx}, ny={args.ny}, nz={args.nz}")
    print(f"Matrix size: {A.shape[0]} x {A.shape[1]}")
    print(f"Symmetry error ||A-A^T||_F: {symmetry_error(A):.3e}")

    res = solve_robin_quader(args.nx, args.ny, args.nz, args.Lx, args.Ly, args.Lz, args.beta, args.modes, args.tol)
    print("\nLowest eigenvalues:")
    for i, (lam, kval) in enumerate(zip(res.eigenvalues, res.k_values), start=1):
        print(f"  mode {i:2d}: lambda = {lam:.12f},   k = {kval:.12f}")

    if args.convergence:
        print("\nConvergence table:")
        table = convergence_table(args.convergence, args.Lx, args.Ly, args.Lz, args.beta, min(args.modes, 8), args.tol)
        for n, vals in table:
            vals_str = ", ".join(f"{v:.8f}" for v in vals)
            print(f"  n={n:3d}: {vals_str}")

    if args.beta_scan_vals:
        print("\nBeta scan:")
        table = beta_scan(args.beta_scan_vals, args.nx, args.ny, args.nz, args.Lx, args.Ly, args.Lz, min(args.modes, 12), args.tol)
        for beta, vals in table:
            vals_str = ", ".join(f"{v:.8f}" for v in vals)
            print(f"  beta={beta:8.4f}: {vals_str}")


if __name__ == "__main__":
    main()
