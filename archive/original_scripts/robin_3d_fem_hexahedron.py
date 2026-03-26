import argparse
import math
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class Mesh:
    points: np.ndarray          # (N,3)
    tets: np.ndarray            # (M,4)
    boundary_tris: np.ndarray   # (K,3)


def trilinear_map(xi, eta, zeta, corners):
    """Map from unit cube [0,1]^3 to hexahedron with 8 corners.
    Corner order:
      0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
      4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
    """
    N = np.array([
        (1-xi)*(1-eta)*(1-zeta),
        xi*(1-eta)*(1-zeta),
        xi*eta*(1-zeta),
        (1-xi)*eta*(1-zeta),
        (1-xi)*(1-eta)*zeta,
        xi*(1-eta)*zeta,
        xi*eta*zeta,
        (1-xi)*eta*zeta,
    ])
    return N @ corners


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

    # 6-tet subdivision of each logical cube along the main body diagonal.
    # local order in cube:
    # v000,v100,v010,v110,v001,v101,v011,v111
    tet_pattern = [
        (0, 1, 3, 7),
        (0, 3, 2, 7),
        (0, 2, 6, 7),
        (0, 6, 4, 7),
        (0, 4, 5, 7),
        (0, 5, 1, 7),
    ]

    tets: List[Tuple[int, int, int, int]] = []

    # Boundary triangles from each boundary logical face, split into 2 tris.
    boundary_faces: List[Tuple[int, int, int]] = []

    def add_face(a, b, c, d):
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

                # Boundary faces
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

    return Mesh(points=points, tets=np.asarray(tets, dtype=int), boundary_tris=np.asarray(boundary_faces, dtype=int))


def tet_element_matrices(coords: np.ndarray):
    # coords shape (4,3)
    A = np.ones((4, 4), dtype=float)
    A[:, 1:] = coords
    detA = np.linalg.det(A)
    volume = abs(detA) / 6.0
    if volume <= 0:
        raise ValueError("Degenerate tetrahedron encountered")

    invA = np.linalg.inv(A)
    # grad phi_i = coefficients of x,y,z in barycentric coordinate lambda_i
    grads = invA[1:, :]  # shape (3,4)
    K = volume * (grads.T @ grads)
    M = (volume / 20.0) * np.array([
        [2, 1, 1, 1],
        [1, 2, 1, 1],
        [1, 1, 2, 1],
        [1, 1, 1, 2],
    ], dtype=float)
    return K, M


def tri_boundary_matrix(coords: np.ndarray):
    # coords shape (3,3)
    area = 0.5 * np.linalg.norm(np.cross(coords[1] - coords[0], coords[2] - coords[0]))
    R = (area / 12.0) * np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ], dtype=float)
    return R


def assemble(mesh: Mesh, beta: float):
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

    A = (K + beta * R).tocsr()
    M = M.tocsr()
    return A, M, R.tocsr()


def solve_modes(A: csr_matrix, M: csr_matrix, modes: int):
    vals, vecs = eigsh(A, M=M, k=modes, sigma=0.0, which='LM')
    order = np.argsort(vals)
    vals = vals[order]
    vecs = vecs[:, order]
    # M-normalize
    for i in range(vecs.shape[1]):
        nrm = math.sqrt(vecs[:, i].T @ (M @ vecs[:, i]))
        vecs[:, i] /= nrm
    return vals, vecs


def make_corners(Lx, Ly, Lz, shear_xy=0.0, shear_xz=0.0, shear_yz=0.0):
    # Start from box and shear top face.
    c = np.array([
        [0.0, 0.0, 0.0],
        [Lx, 0.0, 0.0],
        [Lx, Ly, 0.0],
        [0.0, Ly, 0.0],
        [shear_xz, shear_yz, Lz],
        [Lx + shear_xz, shear_xy + shear_yz, Lz],
        [Lx + shear_xz, Ly + shear_xy + shear_yz, Lz],
        [shear_xz, Ly + shear_yz, Lz],
    ], dtype=float)
    return c


def symmetry_check(A: csr_matrix):
    D = A - A.T
    return np.sqrt((D.multiply(D)).sum())


def run_case(nx, ny, nz, corners, beta, modes):
    mesh = build_hexahedral_mesh(nx, ny, nz, corners)
    A, M, R = assemble(mesh, beta)
    vals, vecs = solve_modes(A, M, modes)
    return mesh, A, M, vals, vecs


def main():
    ap = argparse.ArgumentParser(description='3D Robin FEM eigenproblem on a deformed hexahedron.')
    ap.add_argument('--nx', type=int, default=8)
    ap.add_argument('--ny', type=int, default=8)
    ap.add_argument('--nz', type=int, default=8)
    ap.add_argument('--Lx', type=float, default=1.0)
    ap.add_argument('--Ly', type=float, default=1.0)
    ap.add_argument('--Lz', type=float, default=1.0)
    ap.add_argument('--shear-xy', type=float, default=0.0)
    ap.add_argument('--shear-xz', type=float, default=0.0)
    ap.add_argument('--shear-yz', type=float, default=0.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--modes', type=int, default=8)
    ap.add_argument('--convergence', type=int, nargs='*', default=None,
                    help='Run with isotropic refinements n, using nx=ny=nz=n')
    args = ap.parse_args()

    corners = make_corners(args.Lx, args.Ly, args.Lz, args.shear_xy, args.shear_xz, args.shear_yz)

    if args.convergence:
        print('Convergence study')
        print('n, lambda_1, lambda_2, lambda_3, lambda_4')
        for n in args.convergence:
            _, A, M, vals, _ = run_case(n, n, n, corners, args.beta, min(args.modes, 4))
            print(f'{n}, ' + ', '.join(f'{x:.12f}' for x in vals[:4]))
        return

    mesh, A, M, vals, vecs = run_case(args.nx, args.ny, args.nz, corners, args.beta, args.modes)
    asym = symmetry_check(A)

    print('3D Robin FEM eigenproblem on deformed hexahedron')
    print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
    print(f'beta={args.beta}')
    print(f'symmetry ||A-A^T||_F = {asym:.6e}')
    print('First eigenvalues:')
    for i, lam in enumerate(vals, start=1):
        k = math.sqrt(max(lam, 0.0))
        print(f'  mode {i:2d}: lambda={lam:.12f}, k={k:.12f}')


if __name__ == '__main__':
    main()
