import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class Mesh:
    points: np.ndarray
    tets: np.ndarray
    boundary_tris: np.ndarray


def trilinear_map(xi: float, eta: float, zeta: float, corners: np.ndarray) -> np.ndarray:
    N = np.array([
        (1 - xi) * (1 - eta) * (1 - zeta),
        xi * (1 - eta) * (1 - zeta),
        xi * eta * (1 - zeta),
        (1 - xi) * eta * (1 - zeta),
        (1 - xi) * (1 - eta) * zeta,
        xi * (1 - eta) * zeta,
        xi * eta * zeta,
        (1 - xi) * eta * zeta,
    ])
    return N @ corners


def make_corners(shape: str = "cube", Lx: float = 1.0, Ly: float = 1.0, Lz: float = 1.0,
                 shear_xy: float = 0.0, shear_xz: float = 0.0, shear_yz: float = 0.0) -> np.ndarray:
    shape = shape.lower()
    if shape not in {"cube", "quader", "hexahedron"}:
        raise ValueError(f"unknown shape: {shape}")
    if shape == "cube":
        Lx = Ly = Lz = 1.0
    return np.array([
        [0.0, 0.0, 0.0],
        [Lx, 0.0, 0.0],
        [Lx, Ly, 0.0],
        [0.0, Ly, 0.0],
        [shear_xz, shear_yz, Lz],
        [Lx + shear_xz, shear_xy + shear_yz, Lz],
        [Lx + shear_xz, Ly + shear_xy + shear_yz, Lz],
        [shear_xz, Ly + shear_yz, Lz],
    ], dtype=float)


def build_hexahedral_mesh(nx: int, ny: int, nz: int, corners: np.ndarray) -> Mesh:
    xs = np.linspace(0.0, 1.0, nx)
    ys = np.linspace(0.0, 1.0, ny)
    zs = np.linspace(0.0, 1.0, nz)

    node_id: Dict[Tuple[int, int, int], int] = {}
    points: List[np.ndarray] = []
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                node_id[(i, j, k)] = len(points)
                points.append(trilinear_map(x, y, z, corners))
    points = np.asarray(points, dtype=float)

    tet_pattern = [
        (0, 1, 3, 7),
        (0, 3, 2, 7),
        (0, 2, 6, 7),
        (0, 6, 4, 7),
        (0, 4, 5, 7),
        (0, 5, 1, 7),
    ]
    tets: List[Tuple[int, int, int, int]] = []
    boundary_faces: List[Tuple[int, int, int]] = []

    def add_face(a: int, b: int, c: int, d: int) -> None:
        boundary_faces.append((a, b, c))
        boundary_faces.append((a, c, d))

    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                v000 = node_id[(i, j, k)]
                v100 = node_id[(i + 1, j, k)]
                v010 = node_id[(i, j + 1, k)]
                v110 = node_id[(i + 1, j + 1, k)]
                v001 = node_id[(i, j, k + 1)]
                v101 = node_id[(i + 1, j, k + 1)]
                v011 = node_id[(i, j + 1, k + 1)]
                v111 = node_id[(i + 1, j + 1, k + 1)]
                cube = [v000, v100, v010, v110, v001, v101, v011, v111]
                for tet in tet_pattern:
                    tets.append(tuple(cube[idx] for idx in tet))
                if i == 0:
                    add_face(v000, v010, v011, v001)
                if i == nx - 2:
                    add_face(v100, v101, v111, v110)
                if j == 0:
                    add_face(v000, v001, v101, v100)
                if j == ny - 2:
                    add_face(v010, v110, v111, v011)
                if k == 0:
                    add_face(v000, v100, v110, v010)
                if k == nz - 2:
                    add_face(v001, v011, v111, v101)

    return Mesh(points=points,
                tets=np.asarray(tets, dtype=int),
                boundary_tris=np.asarray(boundary_faces, dtype=int))


def tet_element_matrices(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A = np.ones((4, 4), dtype=float)
    A[:, 1:] = coords
    detA = np.linalg.det(A)
    volume = abs(detA) / 6.0
    if volume <= 0:
        raise ValueError("Degenerate tetrahedron encountered")
    invA = np.linalg.inv(A)
    grads = invA[1:, :]
    K = volume * (grads.T @ grads)
    M = (volume / 20.0) * np.array([
        [2, 1, 1, 1],
        [1, 2, 1, 1],
        [1, 1, 2, 1],
        [1, 1, 1, 2],
    ], dtype=float)
    return K, M


def tri_boundary_matrix(coords: np.ndarray) -> np.ndarray:
    area = 0.5 * np.linalg.norm(np.cross(coords[1] - coords[0], coords[2] - coords[0]))
    return (area / 12.0) * np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ], dtype=float)


def assemble(mesh: Mesh, beta: float) -> tuple[csr_matrix, csr_matrix, csr_matrix]:
    n = mesh.points.shape[0]
    K = lil_matrix((n, n), dtype=float)
    M = lil_matrix((n, n), dtype=float)
    R = lil_matrix((n, n), dtype=float)

    for tet in mesh.tets:
        coords = mesh.points[tet]
        Ke, Me = tet_element_matrices(coords)
        for a in range(4):
            ia = tet[a]
            for b in range(4):
                ib = tet[b]
                K[ia, ib] += Ke[a, b]
                M[ia, ib] += Me[a, b]

    for tri in mesh.boundary_tris:
        coords = mesh.points[tri]
        Re = tri_boundary_matrix(coords)
        for a in range(3):
            ia = tri[a]
            for b in range(3):
                ib = tri[b]
                R[ia, ib] += Re[a, b]

    return (K + beta * R).tocsr(), M.tocsr(), R.tocsr()


def solve_modes(A: csr_matrix, M: csr_matrix, modes: int) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = eigsh(A, M=M, k=modes, sigma=0.0, which='LM')
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    for i in range(vecs.shape[1]):
        nrm = math.sqrt(vecs[:, i].T @ (M @ vecs[:, i]))
        vecs[:, i] /= nrm
    return vals, vecs


def symmetry_norm(A: csr_matrix) -> float:
    D = A - A.T
    return float(np.sqrt((D.multiply(D)).sum()))
