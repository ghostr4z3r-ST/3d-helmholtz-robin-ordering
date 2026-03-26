import argparse
import math
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


@dataclass
class Mesh:
    points: np.ndarray
    tets: np.ndarray
    boundary_tris: np.ndarray
    logical_coords: np.ndarray  # (N,3) logical cube coords in [0,1]^3


def trilinear_map(xi, eta, zeta, corners):
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


def make_corners(kind: str, Lx=1.0, Ly=1.0, Lz=1.0, shear_xy=0.0, shear_xz=0.0, shear_yz=0.0):
    if kind == 'cube':
        Lx = Ly = Lz = 1.0
        shear_xy = shear_xz = shear_yz = 0.0
    elif kind == 'quader':
        shear_xy = shear_xz = shear_yz = 0.0
    elif kind == 'hexahedron':
        pass
    else:
        raise ValueError(f'Unknown geometry kind: {kind}')
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


def build_hexahedral_mesh(nx: int, ny: int, nz: int, corners: np.ndarray) -> Mesh:
    xs = np.linspace(0.0, 1.0, nx)
    ys = np.linspace(0.0, 1.0, ny)
    zs = np.linspace(0.0, 1.0, nz)

    node_id: Dict[Tuple[int, int, int], int] = {}
    points: List[np.ndarray] = []
    logicals: List[np.ndarray] = []

    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                node_id[(i, j, k)] = len(points)
                logicals.append(np.array([x, y, z], dtype=float))
                points.append(trilinear_map(x, y, z, corners))

    points = np.asarray(points, dtype=float)
    logicals = np.asarray(logicals, dtype=float)

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

    return Mesh(points=points, tets=np.asarray(tets, dtype=int), boundary_tris=np.asarray(boundary_faces, dtype=int), logical_coords=logicals)


def tet_element_matrices(coords: np.ndarray):
    A = np.ones((4, 4), dtype=float)
    A[:, 1:] = coords
    detA = np.linalg.det(A)
    volume = abs(detA) / 6.0
    if volume <= 0:
        raise ValueError('Degenerate tetrahedron encountered')
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


def tri_boundary_matrix(coords: np.ndarray):
    area = 0.5 * np.linalg.norm(np.cross(coords[1] - coords[0], coords[2] - coords[0]))
    return (area / 12.0) * np.array([
        [2, 1, 1],
        [1, 2, 1],
        [1, 1, 2],
    ], dtype=float)


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
    return A, M.tocsr()


def solve_modes(A: csr_matrix, M: csr_matrix, modes: int):
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
    return math.sqrt((D.multiply(D)).sum())


# --- Local 8-octant diagnostics -------------------------------------------------

def hadamard_basis_8():
    sigs = np.array([
        [-1, -1, -1], [ 1, -1, -1], [-1,  1, -1], [ 1,  1, -1],
        [-1, -1,  1], [ 1, -1,  1], [-1,  1,  1], [ 1,  1,  1],
    ], dtype=int)
    names = ['const', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']
    basis = []
    for name in names:
        if name == 'const':
            v = np.ones(8)
        elif name == 'x':
            v = sigs[:, 0]
        elif name == 'y':
            v = sigs[:, 1]
        elif name == 'z':
            v = sigs[:, 2]
        elif name == 'xy':
            v = sigs[:, 0] * sigs[:, 1]
        elif name == 'xz':
            v = sigs[:, 0] * sigs[:, 2]
        elif name == 'yz':
            v = sigs[:, 1] * sigs[:, 2]
        elif name == 'xyz':
            v = sigs[:, 0] * sigs[:, 1] * sigs[:, 2]
        v = v.astype(float)
        v /= np.linalg.norm(v)
        basis.append(v)
    return names, np.vstack(basis), sigs


def local_octant_signature(mesh: Mesh, mode: np.ndarray, logical_center: Tuple[float, float, float], radius=0.22):
    lc = np.asarray(logical_center, dtype=float)
    dx = mesh.logical_coords - lc[None, :]
    mask = np.all(np.abs(dx) <= radius, axis=1)
    # exclude points too close to splitting planes to avoid ambiguity
    eps = 1e-12
    mask &= np.all(np.abs(dx) > eps, axis=1)
    selected = np.where(mask)[0]
    vals8 = np.zeros(8, dtype=float)
    counts = np.zeros(8, dtype=int)
    names, basis, sigs = hadamard_basis_8()
    if len(selected) == 0:
        raise RuntimeError('No local points selected; increase radius or resolution.')
    for idx in selected:
        s = np.sign(dx[idx]).astype(int)
        # map sign triplet to basis ordering above
        oct_idx = np.where((sigs == s).all(axis=1))[0][0]
        vals8[oct_idx] += mode[idx]
        counts[oct_idx] += 1
    for i in range(8):
        if counts[i] > 0:
            vals8[i] /= counts[i]
    # if some octants empty, fail loudly
    if np.any(counts == 0):
        raise RuntimeError(f'Local diagnostic failed: some octants empty at {logical_center}, counts={counts.tolist()}')
    vals8 -= vals8.mean()  # remove constant offset to emphasize inner sign structure
    nrm = np.linalg.norm(vals8)
    if nrm > 0:
        vals8n = vals8 / nrm
    else:
        vals8n = vals8
    coeffs = basis @ vals8n
    order = np.argsort(-np.abs(coeffs))
    ranked = [(names[i], float(coeffs[i])) for i in order]
    return {
        'values8': vals8,
        'counts': counts,
        'coeffs': dict((names[i], float(coeffs[i])) for i in range(8)),
        'ranked': ranked,
    }


def first_xyz_mode(mesh: Mesh, vecs: np.ndarray, logical_center=(0.5, 0.5, 0.5), radius=0.22):
    best = None
    for i in range(vecs.shape[1]):
        sig = local_octant_signature(mesh, vecs[:, i], logical_center, radius)
        xyz = abs(sig['coeffs']['xyz'])
        if best is None or xyz > best[1]:
            best = (i, xyz, sig)
    return best


def print_mode_summary(vals, vecs, mesh, mode_indices, locations, radius=0.22):
    for mi in mode_indices:
        idx = mi - 1
        print(f"Mode {mi}: lambda={vals[idx]:.12f}")
        for loc_name, loc in locations.items():
            sig = local_octant_signature(mesh, vecs[:, idx], loc, radius=radius)
            top3 = ', '.join([f"{n}:{c:+.3f}" for n, c in sig['ranked'][:3]])
            print(f"  {loc_name:<10s} top signatures -> {top3}")


def main():
    ap = argparse.ArgumentParser(description='Center-emergence diagnostics for 3D Robin FEM on cube/quader/hexahedron.')
    ap.add_argument('--geometry', choices=['cube', 'quader', 'hexahedron'], default='cube')
    ap.add_argument('--nx', type=int, default=10)
    ap.add_argument('--ny', type=int, default=10)
    ap.add_argument('--nz', type=int, default=10)
    ap.add_argument('--Lx', type=float, default=1.0)
    ap.add_argument('--Ly', type=float, default=1.0)
    ap.add_argument('--Lz', type=float, default=1.0)
    ap.add_argument('--shear-xy', type=float, default=0.0)
    ap.add_argument('--shear-xz', type=float, default=0.0)
    ap.add_argument('--shear-yz', type=float, default=0.0)
    ap.add_argument('--beta', type=float, default=1.0)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--radius', type=float, default=0.22, help='Local box half-width in logical coordinates.')
    ap.add_argument('--mode-summary', type=int, nargs='*', default=None,
                    help='Mode numbers for detailed local signature summaries.')
    args = ap.parse_args()

    corners = make_corners(args.geometry, args.Lx, args.Ly, args.Lz, args.shear_xy, args.shear_xz, args.shear_yz)
    mesh = build_hexahedral_mesh(args.nx, args.ny, args.nz, corners)
    A, M = assemble(mesh, args.beta)
    vals, vecs = solve_modes(A, M, args.modes)

    print(f'Geometry={args.geometry}, beta={args.beta}')
    print(f'nodes={mesh.points.shape[0]}, tets={mesh.tets.shape[0]}, boundary_tris={mesh.boundary_tris.shape[0]}')
    print(f'symmetry ||A-A^T||_F = {symmetry_norm(A):.6e}')
    print('First eigenvalues:')
    for i, lam in enumerate(vals, start=1):
        print(f'  mode {i:2d}: lambda={lam:.12f}')

    locations = {
        'center': (0.5, 0.5, 0.5),
        'face_x0': (0.12, 0.5, 0.5),
        'edge_x0y0': (0.12, 0.12, 0.5),
        'corner_near': (0.12, 0.12, 0.12),
    }

    mode_idx, xyz_amp, sig = first_xyz_mode(mesh, vecs, logical_center=locations['center'], radius=args.radius)
    print('\nEarliest/strongest local xyz-like center signature among computed modes:')
    print(f'  best mode = {mode_idx+1}, lambda={vals[mode_idx]:.12f}, |xyz coeff|={xyz_amp:.6f}')
    print('  center top signatures: ' + ', '.join([f"{n}:{c:+.3f}" for n, c in sig['ranked'][:5]]))

    if args.mode_summary:
        print('\nDetailed local signature summaries')
        print_mode_summary(vals, vecs, mesh, args.mode_summary, locations, radius=args.radius)


if __name__ == '__main__':
    main()
