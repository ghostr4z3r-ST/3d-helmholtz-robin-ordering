from __future__ import annotations

from typing import Dict, Iterable, Tuple

import numpy as np


def nearest_indices(points: np.ndarray, center: np.ndarray, count: int = 64) -> np.ndarray:
    d = np.linalg.norm(points - center[None, :], axis=1)
    return np.argsort(d)[: min(count, len(points))]


def bbox_normalize(points: np.ndarray, query: Iterable[float]) -> np.ndarray:
    p = np.asarray(query, dtype=float)
    mins = points.min(axis=0)
    maxs = points.max(axis=0)
    span = np.where(maxs > mins, maxs - mins, 1.0)
    return mins + p * span


def local_octant_signature(mesh, vec: np.ndarray, logical_center: Tuple[float, float, float] | None = None,
                           physical_center: Tuple[float, float, float] | None = None,
                           radius: float = 0.22, count: int = 128) -> Dict[str, object]:
    if (logical_center is None) == (physical_center is None):
        raise ValueError("provide exactly one of logical_center or physical_center")

    if logical_center is not None:
        center = bbox_normalize(mesh.points, logical_center)
        span = mesh.points.max(axis=0) - mesh.points.min(axis=0)
        scale = float(np.max(span)) if np.max(span) > 0 else 1.0
        radius_phys = radius * scale
    else:
        center = np.asarray(physical_center, dtype=float)
        radius_phys = radius

    d = np.linalg.norm(mesh.points - center[None, :], axis=1)
    idx = np.where(d <= radius_phys)[0]
    if idx.size < 8:
        idx = nearest_indices(mesh.points, center, count=count)

    rel = mesh.points[idx] - center[None, :]
    vals = vec[idx].astype(float).copy()
    vals -= vals.mean()

    patterns = [(-1, -1, -1), (1, -1, -1), (-1, 1, -1), (1, 1, -1),
                (-1, -1, 1), (1, -1, 1), (-1, 1, 1), (1, 1, 1)]
    oct_means = []
    for sx, sy, sz in patterns:
        mask = ((np.sign(rel[:, 0] + 1e-15) == sx)
                & (np.sign(rel[:, 1] + 1e-15) == sy)
                & (np.sign(rel[:, 2] + 1e-15) == sz))
        if np.any(mask):
            oct_means.append(float(vals[mask].mean()))
        else:
            dirs = np.array([sx, sy, sz], dtype=float)
            score = rel @ dirs
            j = int(np.argmax(score))
            oct_means.append(float(vals[j]))
    q = np.asarray(oct_means, dtype=float)

    basis = {
        'const': np.array([1, 1, 1, 1, 1, 1, 1, 1], float),
        'x':     np.array([-1, 1, -1, 1, -1, 1, -1, 1], float),
        'y':     np.array([-1, -1, 1, 1, -1, -1, 1, 1], float),
        'z':     np.array([-1, -1, -1, -1, 1, 1, 1, 1], float),
        'xy':    np.array([1, -1, -1, 1, 1, -1, -1, 1], float),
        'xz':    np.array([1, -1, 1, -1, -1, 1, -1, 1], float),
        'yz':    np.array([1, 1, -1, -1, -1, -1, 1, 1], float),
        'xyz':   np.array([-1, 1, 1, -1, 1, -1, -1, 1], float),
    }

    qn = float(np.linalg.norm(q))
    coeffs: Dict[str, float] = {}
    for name, b in basis.items():
        bn = float(np.linalg.norm(b))
        coeffs[name] = float(np.dot(q, b) / (qn * bn)) if qn > 0 else 0.0

    ranked = sorted(coeffs.items(), key=lambda kv: abs(kv[1]), reverse=True)
    return {
        'center': center,
        'radius_phys': radius_phys,
        'octant_means': q,
        'coeffs': coeffs,
        'ranked': ranked,
        'support_size': int(idx.size),
    }


def first_xyz_mode(mesh, vecs: np.ndarray, logical_center=(0.5, 0.5, 0.5), radius: float = 0.22):
    best_idx = None
    best_val = -1.0
    best_sig = None
    for i in range(vecs.shape[1]):
        sig = local_octant_signature(mesh, vecs[:, i], logical_center=logical_center, radius=radius)
        val = abs(sig['coeffs']['xyz'])
        if val > best_val:
            best_idx = i
            best_val = val
            best_sig = sig
    return int(best_idx), float(best_val), best_sig
