
"""
3D Robin eigenvalue problem on a cube using tensor-product linear finite elements.

We solve
    -Δu = λ u  in Ω = [0, L]^3
with Robin boundary condition
    ∂_n u + beta * u = 0  on ∂Ω.

Weak form:
    ∫_Ω ∇u·∇v dΩ + beta ∫_{∂Ω} u v dS = λ ∫_Ω u v dΩ

This implementation builds the symmetric generalized eigenproblem
    A u = λ M u
using tensor-product 1D FEM matrices on a uniform grid.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class EigenResult:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    k_values: np.ndarray


def fem1d_matrices(n: int, L: float = 1.0) -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix]:
    """
    1D linear FEM matrices on [0, L] with n nodes.
    Returns:
        K: stiffness matrix for ∫ u' v'
        M: mass matrix for ∫ u v
        B: boundary matrix for endpoint contributions v(0)u(0) + v(L)u(L)
    """
    if n < 2:
        raise ValueError("n must be at least 2.")

    h = L / (n - 1)

    # Stiffness matrix
    main_k = np.full(n, 2.0 / h)
    off_k = np.full(n - 1, -1.0 / h)
    K = sp.diags([off_k, main_k, off_k], offsets=[-1, 0, 1], format="lil")
    K[0, 0] = 1.0 / h
    K[-1, -1] = 1.0 / h
    K = K.tocsr()

    # Mass matrix
    main_m = np.full(n, 4.0 * h / 6.0)
    off_m = np.full(n - 1, 1.0 * h / 6.0)
    M = sp.diags([off_m, main_m, off_m], offsets=[-1, 0, 1], format="lil")
    M[0, 0] = h / 3.0
    M[-1, -1] = h / 3.0
    M = M.tocsr()

    # Boundary endpoint matrix
    bdiag = np.zeros(n, dtype=float)
    bdiag[0] = 1.0
    bdiag[-1] = 1.0
    B = sp.diags(bdiag, 0, format="csr")

    return K, M, B


def robin_3d_matrices(n: int, L: float = 1.0, beta: float = 1.0) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """
    Build the 3D FEM matrices on the cube [0, L]^3 using tensor products.
    Returns:
        A: stiffness + Robin boundary matrix
        M: mass matrix
    """
    K1, M1, B1 = fem1d_matrices(n=n, L=L)

    # Volume terms: ∫ (ux vx + uy vy + uz vz)
    A_vol = (
        sp.kron(sp.kron(K1, M1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, K1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, M1, format="csr"), K1, format="csr")
    )

    # Robin boundary terms on the six cube faces
    R = beta * (
        sp.kron(sp.kron(B1, M1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, B1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, M1, format="csr"), B1, format="csr")
    )

    M = sp.kron(sp.kron(M1, M1, format="csr"), M1, format="csr")

    return (A_vol + R).tocsr(), M.tocsr()


def solve_robin_cube(
    n: int = 10,
    L: float = 1.0,
    beta: float = 1.0,
    modes: int = 10,
    tol: float = 1e-9,
) -> EigenResult:
    """
    Solve the generalized symmetric eigenproblem A u = λ M u for the lowest modes.
    """
    A, M = robin_3d_matrices(n=n, L=L, beta=beta)
    if modes >= A.shape[0] - 1:
        raise ValueError("modes must be smaller than the total number of degrees of freedom minus 1.")

    vals, vecs = spla.eigsh(A, k=modes, M=M, sigma=0.0, which="LM", tol=tol)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    kvals = np.sqrt(np.clip(vals, 0.0, None))
    return EigenResult(eigenvalues=vals, eigenvectors=vecs, k_values=kvals)


def symmetry_error(A: sp.csr_matrix) -> float:
    D = A - A.T
    if D.nnz == 0:
        return 0.0
    return float(np.sqrt((D.data ** 2).sum()))


def convergence_table(ns: list[int], L: float, beta: float, modes: int, tol: float) -> list[tuple[int, np.ndarray]]:
    table = []
    for n in ns:
        res = solve_robin_cube(n=n, L=L, beta=beta, modes=modes, tol=tol)
        table.append((n, res.eigenvalues.copy()))
    return table


def main() -> None:
    parser = argparse.ArgumentParser(description="Solve the 3D Robin eigenvalue problem on a cube using tensor-product FEM.")
    parser.add_argument("--n", type=int, default=10, help="Number of nodes per axis.")
    parser.add_argument("--L", type=float, default=1.0, help="Cube side length.")
    parser.add_argument("--beta", type=float, default=1.0, help="Robin boundary parameter.")
    parser.add_argument("--modes", type=int, default=10, help="Number of smallest eigenmodes to compute.")
    parser.add_argument("--tol", type=float, default=1e-9, help="Eigensolver tolerance.")
    parser.add_argument("--convergence", nargs="*", type=int, help="Optional list of n values for a convergence table.")
    args = parser.parse_args()

    A, M = robin_3d_matrices(n=args.n, L=args.L, beta=args.beta)
    print(f"Matrix size: {A.shape[0]} x {A.shape[1]}")
    print(f"Symmetry error ||A-A^T||_F: {symmetry_error(A):.3e}")

    res = solve_robin_cube(n=args.n, L=args.L, beta=args.beta, modes=args.modes, tol=args.tol)

    print("\nLowest eigenvalues:")
    for i, (lam, kval) in enumerate(zip(res.eigenvalues, res.k_values), start=1):
        print(f"  mode {i:2d}: lambda = {lam:.12f},   k = {kval:.12f}")

    if args.convergence:
        print("\nConvergence table:")
        table = convergence_table(args.convergence, L=args.L, beta=args.beta, modes=min(args.modes, 6), tol=args.tol)
        for n, vals in table:
            vals_str = ", ".join(f"{v:.8f}" for v in vals[: min(6, len(vals))])
            print(f"  n={n:3d}: {vals_str}")


if __name__ == "__main__":
    main()
