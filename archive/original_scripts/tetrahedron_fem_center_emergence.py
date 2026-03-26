import argparse, math
import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh


def regular_tetra_vertices(scale=1.0):
    v = np.array([
        [1.0, 1.0, 1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [1.0, -1.0, -1.0],
    ])
    return scale * v / math.sqrt(3.0)


def barycentric_lattice_points(vertices, n):
    pts = []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                l = n - i - j - k
                w = np.array([i, j, k, l], dtype=float) / n
                p = w @ vertices
                pts.append(p)
    pts = np.unique(np.round(np.array(pts), 12), axis=0)
    return pts


def tetra_volume(coords):
    a, b, c, d = coords
    return abs(np.linalg.det(np.column_stack((b - a, c - a, d - a)))) / 6.0


def tet_gradients(coords):
    T = np.ones((4, 4))
    T[:, 1:] = coords
    invT = np.linalg.inv(T)
    return invT[1:, :]  # columns are grad phi_i


def local_tet_matrices(coords):
    V = tetra_volume(coords)
    grads = tet_gradients(coords)
    K = V * (grads.T @ grads)
    M = (V / 20.0) * (np.ones((4, 4)) + np.eye(4))
    return K, M


def local_tri_mass(coords3):
    a, b, c = coords3
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    return (area / 12.0) * (np.ones((3, 3)) + np.eye(3))


def assemble(vertices, n, beta):
    pts = barycentric_lattice_points(vertices, n)
    dela = Delaunay(pts)
    simplices = dela.simplices.copy()

    N = len(pts)
    K = lil_matrix((N, N))
    M = lil_matrix((N, N))
    R = lil_matrix((N, N))

    # volume terms
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

    # boundary faces via convex hull triangles
    for tri in dela.convex_hull:
        coords3 = pts[tri]
        Re = local_tri_mass(coords3)
        for a in range(3):
            ia = tri[a]
            for b in range(3):
                ib = tri[b]
                R[ia, ib] += Re[a, b]

    A = (K + beta * R).tocsr()
    B = M.tocsr()
    return pts, simplices, A, B


def solve(vertices, n, beta, modes):
    pts, simplices, A, B = assemble(vertices, n, beta)
    vals, vecs = eigsh(A, k=modes, M=B, sigma=0.0, which='LM')
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    return pts, simplices, A, B, vals, vecs


def nearest_indices(points, center, count=27):
    d = np.linalg.norm(points - center, axis=1)
    return np.argsort(d)[:count]


def local_octant_signature(points, vec, center, count=64):
    idx = nearest_indices(points, center, count=min(count, len(points)))
    rel = points[idx] - center
    vals = vec[idx].copy()
    vals = vals - vals.mean()
    # octant averages with fallback nearest point weighting
    oct_means = []
    patterns = [(-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1), (-1,-1,1), (1,-1,1), (-1,1,1), (1,1,1)]
    for sx, sy, sz in patterns:
        mask = (np.sign(rel[:,0] + 1e-15) == sx) & (np.sign(rel[:,1] + 1e-15) == sy) & (np.sign(rel[:,2] + 1e-15) == sz)
        if np.any(mask):
            oct_means.append(vals[mask].mean())
        else:
            # nearest point in that direction by directional score
            dirs = np.array([sx, sy, sz])
            score = rel @ dirs
            j = np.argmax(score)
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


def best_xyz_mode(points, vecs, center):
    best = None
    bestv = -1
    for m in range(vecs.shape[1]):
        coeffs = local_octant_signature(points, vecs[:, m], center)
        v = abs(coeffs['xyz'])
        if v > bestv:
            bestv = v
            best = (m + 1, coeffs)
    return best, bestv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=10)
    ap.add_argument('--beta', type=float, default=0.0)
    ap.add_argument('--modes', type=int, default=12)
    ap.add_argument('--scale', type=float, default=1.0)
    args = ap.parse_args()

    verts = regular_tetra_vertices(args.scale)
    center = verts.mean(axis=0)
    points, simplices, A, B, vals, vecs = solve(verts, args.n, args.beta, args.modes)
    asym = (A - A.T).power(2).sum() ** 0.5
    print(f'asymmetry_norm={asym:.3e}')
    for i, lam in enumerate(vals, start=1):
        print(f'mode {i:2d}: lambda={lam:.12f}')
    (best_mode, coeffs), bestv = best_xyz_mode(points, vecs, center)
    print(f'best_xyz_mode={best_mode}, abs_c_xyz={bestv:.6f}')
    print('center_signature=' + ', '.join(f'{k}:{coeffs[k]:+.3f}' for k in ['const','x','y','z','xy','xz','yz','xyz']))

if __name__ == '__main__':
    main()
