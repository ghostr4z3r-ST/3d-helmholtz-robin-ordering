from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class TetraResult:
    points: np.ndarray
    simplices: np.ndarray
    A: csr_matrix
    M: csr_matrix
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray


def regular_tetra_vertices(scale: float = 1.0) -> np.ndarray:
    v = np.array([
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ])
    return scale * v / math.sqrt(3.0)


def barycentric_lattice_points(vertices: np.ndarray, n: int) -> np.ndarray:
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                l = n - i - j - k
                w = np.array([i, j, k, l], dtype=float) / n
                pts.append(w @ vertices)
    pts = np.unique(np.round(np.array(pts), 12), axis=0)
    return pts


def tetra_volume(coords: np.ndarray) -> float:
    a, b, c, d = coords
    return abs(np.linalg.det(np.column_stack((b - a, c - a, d - a)))) / 6.0


def tet_gradients(coords: np.ndarray) -> np.ndarray:
    T = np.ones((4, 4))
    T[:, 1:] = coords
    invT = np.linalg.inv(T)
    return invT[1:, :]


def local_tet_matrices(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    V = tetra_volume(coords)
    grads = tet_gradients(coords)
    K = V * (grads.T @ grads)
    M = (V / 20.0) * (np.ones((4, 4)) + np.eye(4))
    return K, M


def local_tri_mass(coords3: np.ndarray) -> np.ndarray:
    a, b, c = coords3
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    return (area / 12.0) * (np.ones((3, 3)) + np.eye(3))


def assemble(vertices: np.ndarray, n: int, beta: float) -> tuple[np.ndarray, np.ndarray, csr_matrix, csr_matrix]:
    pts = barycentric_lattice_points(vertices, n)
    dela = Delaunay(pts)
    simplices = dela.simplices.copy()

    N = len(pts)
    K = lil_matrix((N, N))
    M = lil_matrix((N, N))
    R = lil_matrix((N, N))

    for tet in simplices:
        coords = pts[tet]
        V = tetra_volume(coords)
        if V < 1e-12:
            continue
        Ke, Me = local_tet_matrices(coords)
        for a in range(4):
            ia = tet[a]
            for b in range(4):
                ib = tet[b]
                K[ia, ib] += Ke[a, b]
                M[ia, ib] += Me[a, b]

    for tri in dela.convex_hull:
        coords3 = pts[tri]
        Re = local_tri_mass(coords3)
        for a in range(3):
            ia = tri[a]
            for b in range(3):
                ib = tri[b]
                R[ia, ib] += Re[a, b]

    return pts, simplices, (K + beta * R).tocsr(), M.tocsr()


def solve(vertices: np.ndarray, n: int, beta: float, modes: int) -> TetraResult:
    pts, simplices, A, M = assemble(vertices, n, beta)
    vals, vecs = eigsh(A, k=modes, M=M, sigma=0.0, which='LM')
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    return TetraResult(pts, simplices, A, M, vals, vecs)


def local_octant_signature(points: np.ndarray, vec: np.ndarray, center: np.ndarray, count: int = 64) -> dict[str, float]:
    d = np.linalg.norm(points - center, axis=1)
    idx = np.argsort(d)[: min(count, len(points))]
    rel = points[idx] - center
    vals = vec[idx].copy() - vec[idx].mean()
    patterns = [(-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1), (-1,-1,1), (1,-1,1), (-1,1,1), (1,1,1)]
    oct_means = []
    for sx, sy, sz in patterns:
        mask = (np.sign(rel[:,0] + 1e-15) == sx) & (np.sign(rel[:,1] + 1e-15) == sy) & (np.sign(rel[:,2] + 1e-15) == sz)
        if np.any(mask):
            oct_means.append(vals[mask].mean())
        else:
            dirs = np.array([sx, sy, sz])
            j = np.argmax(rel @ dirs)
            oct_means.append(vals[j])
    q = np.array(oct_means)
    basis = {
        'const': np.array([1,1,1,1,1,1,1,1],float),
        'x':     np.array([-1,1,-1,1,-1,1,-1,1],float),
        'y':     np.array([-1,-1,1,1,-1,-1,1,1],float),
        'z':     np.array([-1,-1,-1,-1,1,1,1,1],float),
        'xy':    np.array([1,-1,-1,1,1,-1,-1,1],float),
        'xz':    np.array([1,-1,1,-1,-1,1,-1,1],float),
        'yz':    np.array([1,1,-1,-1,-1,-1,1,1],float),
        'xyz':   np.array([-1,1,1,-1,1,-1,-1,1],float),
    }
    qn = np.linalg.norm(q)
    coeffs = {}
    for name, b in basis.items():
        bn = np.linalg.norm(b)
        coeffs[name] = float(np.dot(q, b) / (qn * bn)) if qn > 0 else 0.0
    return coeffs
