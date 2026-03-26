from __future__ import annotations

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


def fem1d_matrices(n: int, L: float = 1.0) -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix]:
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


def robin_3d_matrices(n: int, L: float = 1.0, beta: float = 1.0) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    K1, M1, B1 = fem1d_matrices(n=n, L=L)
    A_vol = (
        sp.kron(sp.kron(K1, M1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, K1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, M1, format="csr"), K1, format="csr")
    )
    R = beta * (
        sp.kron(sp.kron(B1, M1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, B1, format="csr"), M1, format="csr")
        + sp.kron(sp.kron(M1, M1, format="csr"), B1, format="csr")
    )
    M = sp.kron(sp.kron(M1, M1, format="csr"), M1, format="csr")
    return (A_vol + R).tocsr(), M.tocsr()


def solve_robin_cube(n: int = 10, L: float = 1.0, beta: float = 1.0, modes: int = 10, tol: float = 1e-9) -> EigenResult:
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


def convergence_table(ns: Iterable[int], L: float, beta: float, modes: int, tol: float) -> list[tuple[int, np.ndarray]]:
    table = []
    for n in ns:
        res = solve_robin_cube(n=n, L=L, beta=beta, modes=modes, tol=tol)
        table.append((n, res.eigenvalues.copy()))
    return table


def beta_scan(betas: Iterable[float], n: int, L: float, modes: int, tol: float) -> list[tuple[float, np.ndarray]]:
    rows = []
    for beta in betas:
        res = solve_robin_cube(n=n, L=L, beta=beta, modes=modes, tol=tol)
        rows.append((float(beta), res.eigenvalues.copy()))
    return rows
