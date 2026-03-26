from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# --- Minimal tensor-product FEM cube solver ------------------------------------

def fem1d_matrices(n: int, L: float = 1.0) -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix]:
    if n < 2:
        raise ValueError("n must be at least 2")
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


def solve_robin_cube(n: int, beta: float, modes: int, tol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    A, M = robin_3d_matrices(n=n, L=1.0, beta=beta)
    vals, vecs = spla.eigsh(A, k=modes, M=M, sigma=0.0, which="LM", tol=tol)
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    # M-normalize
    for j in range(vecs.shape[1]):
        nrm = np.sqrt(vecs[:, j].T @ (M @ vecs[:, j]))
        vecs[:, j] /= nrm
    return vals, vecs


# --- Diagnostics ----------------------------------------------------------------

@dataclass
class ModeMetrics:
    beta: float
    family: str
    mode_index: int
    eigenvalue: float
    feature_label: str
    feature_strength: float
    f_eff: float
    entropy_norm: float
    anisotropy_ratio: float
    eig1: float
    eig2: float
    eig3: float
    boundary_ratio: float


def trap_weights_1d(n: int, L: float = 1.0) -> np.ndarray:
    h = L / (n - 1)
    w = np.ones(n)
    w[0] = 0.5
    w[-1] = 0.5
    return w * h


def make_basis(n: int) -> Tuple[np.ndarray, List[str]]:
    x = np.linspace(0.0, 1.0, n)
    xc = x - 0.5
    X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")
    raws = {
        "const": np.ones_like(X),
        "x": X,
        "y": Y,
        "z": Z,
        "xy": X * Y,
        "xz": X * Z,
        "yz": Y * Z,
        "xyz": X * Y * Z,
    }
    w1 = trap_weights_1d(n)
    W = w1[:, None, None] * w1[None, :, None] * w1[None, None, :]
    basis = []
    names = []
    for name, arr in raws.items():
        nrm = np.sqrt(np.sum(W * arr * arr))
        arrn = arr / nrm
        basis.append(arrn)
        names.append(name)
    return np.stack(basis, axis=0), names


def classify_modes(vals: np.ndarray, vecs: np.ndarray, n: int, max_considered: int = 12) -> Dict[str, Tuple[int, str, float]]:
    basis, names = make_basis(n)
    w1 = trap_weights_1d(n)
    W = w1[:, None, None] * w1[None, :, None] * w1[None, None, :]
    best = {"axis": (-1, "", -1.0), "pair": (-1, "", -1.0), "xyz": (-1, "", -1.0)}
    for j in range(min(max_considered, vecs.shape[1])):
        u = vecs[:, j].reshape((n, n, n))
        # normalize in trapezoidal metric for diagnostics
        u = u / np.sqrt(np.sum(W * u * u))
        coeffs = np.array([abs(np.sum(W * u * b)) for b in basis])
        family_groups = {
            "axis": [1, 2, 3],
            "pair": [4, 5, 6],
            "xyz": [7],
        }
        for fam, idxs in family_groups.items():
            local = coeffs[idxs]
            k = idxs[int(np.argmax(local))]
            strength = float(coeffs[k])
            if strength > best[fam][2]:
                best[fam] = (j + 1, names[k], strength)
    return best


def mode_distribution_metrics(u: np.ndarray, beta: float, family: str, mode_index: int, feature_label: str, feature_strength: float) -> ModeMetrics:
    n = u.shape[0]
    x = np.linspace(0.0, 1.0, n)
    w1 = trap_weights_1d(n)
    W = w1[:, None, None] * w1[None, :, None] * w1[None, None, :]

    # L2-normalize in trapezoidal metric
    u = u / np.sqrt(np.sum(W * u * u))
    rho = u * u
    rho = rho / np.sum(W * rho)

    ipr = float(np.sum(W * rho * rho))
    f_eff = float(1.0 / ipr)  # domain volume is 1

    # Discrete Shannon entropy with normalized voxel probabilities
    p = (W * rho).ravel()
    p = p / p.sum()
    eps = 1e-300
    H = float(-np.sum(p * np.log(p + eps)))
    H_norm = float(H / np.log(len(p)))

    # Second moment tensor
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    xbar = float(np.sum(W * rho * X))
    ybar = float(np.sum(W * rho * Y))
    zbar = float(np.sum(W * rho * Z))
    DX = X - xbar
    DY = Y - ybar
    DZ = Z - zbar
    M = np.array([
        [np.sum(W * rho * DX * DX), np.sum(W * rho * DX * DY), np.sum(W * rho * DX * DZ)],
        [np.sum(W * rho * DY * DX), np.sum(W * rho * DY * DY), np.sum(W * rho * DY * DZ)],
        [np.sum(W * rho * DZ * DX), np.sum(W * rho * DZ * DY), np.sum(W * rho * DZ * DZ)],
    ], dtype=float)
    eigs = np.sort(np.linalg.eigvalsh(M))[::-1]
    anis = float(eigs[-1] / eigs[0]) if eigs[0] > 0 else 0.0

    # Boundary ratio using |u|^2 surface integrals on 6 faces
    wx = w1[:, None]
    # yz faces x=0,1
    face_x0 = np.sum(wx * rho[0, :, :])
    face_x1 = np.sum(wx * rho[-1, :, :])
    # xz faces y=0,1
    face_y0 = np.sum(wx * rho[:, 0, :])
    face_y1 = np.sum(wx * rho[:, -1, :])
    # xy faces z=0,1
    face_z0 = np.sum(wx * rho[:, :, 0])
    face_z1 = np.sum(wx * rho[:, :, -1])
    boundary_ratio = float(face_x0 + face_x1 + face_y0 + face_y1 + face_z0 + face_z1)

    return ModeMetrics(
        beta=beta,
        family=family,
        mode_index=mode_index,
        eigenvalue=np.nan,
        feature_label=feature_label,
        feature_strength=feature_strength,
        f_eff=f_eff,
        entropy_norm=H_norm,
        anisotropy_ratio=anis,
        eig1=float(eigs[0]),
        eig2=float(eigs[1]),
        eig3=float(eigs[2]),
        boundary_ratio=boundary_ratio,
    )


def run_first_test(betas: List[float], n: int, modes: int) -> List[ModeMetrics]:
    rows: List[ModeMetrics] = []
    for beta in betas:
        vals, vecs = solve_robin_cube(n=n, beta=beta, modes=modes)
        best = classify_modes(vals, vecs, n=n, max_considered=min(12, modes))
        for family in ["axis", "pair", "xyz"]:
            mode_index, feature_label, feature_strength = best[family]
            u = vecs[:, mode_index - 1].reshape((n, n, n))
            metrics = mode_distribution_metrics(u, beta, family, mode_index, feature_label, feature_strength)
            metrics.eigenvalue = float(vals[mode_index - 1])
            rows.append(metrics)
    return rows


def write_csv(rows: List[ModeMetrics], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "beta", "family", "mode_index", "eigenvalue", "feature_label", "feature_strength",
            "f_eff", "entropy_norm", "anisotropy_ratio", "eig1", "eig2", "eig3", "boundary_ratio"
        ])
        for r in rows:
            writer.writerow([
                r.beta, r.family, r.mode_index, r.eigenvalue, r.feature_label, r.feature_strength,
                r.f_eff, r.entropy_norm, r.anisotropy_ratio, r.eig1, r.eig2, r.eig3, r.boundary_ratio
            ])


def summarize(rows: List[ModeMetrics]) -> str:
    lines = []
    for beta in sorted({r.beta for r in rows}):
        lines.append(f"beta={beta}")
        sub = [r for r in rows if r.beta == beta]
        for fam in ["axis", "pair", "xyz"]:
            r = next(rr for rr in sub if rr.family == fam)
            lines.append(
                f"  {fam:4s}: mode {r.mode_index:2d}, feat={r.feature_label:>3s} ({r.feature_strength:.3f}), "
                f"f_eff={r.f_eff:.4f}, H={r.entropy_norm:.4f}, anis={r.anisotropy_ratio:.4f}, B={r.boundary_ratio:.4f}"
            )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="First test for field-information distribution in 3D Helmholtz-Robin cube modes.")
    parser.add_argument("--n", type=int, default=17)
    parser.add_argument("--modes", type=int, default=12)
    parser.add_argument("--betas", type=str, default="0,1,5")
    parser.add_argument("--csv", type=str, default="/mnt/data/field_information_distribution_first_test.csv")
    args = parser.parse_args()
    betas = [float(x) for x in args.betas.split(",") if x.strip()]
    rows = run_first_test(betas, n=args.n, modes=args.modes)
    out = Path(args.csv)
    write_csv(rows, out)
    print(summarize(rows))
    print(f"\nSaved CSV to: {out}")


if __name__ == "__main__":
    main()
