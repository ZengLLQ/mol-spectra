"""
Chemical similarity and distance metrics for molecular fingerprints.
"""

import math
from typing import Dict, List, Optional, Tuple
from .fingerprint import MorganFingerprint


def _get_intersection_counts(fp1: MorganFingerprint, fp2: MorganFingerprint) -> Tuple[int, int, int]:
    """
    Returns (a, b, c) where:
    a = popcount(fp1)
    b = popcount(fp2)
    c = popcount(fp1 & fp2)
    """
    if fp1.n_bits != fp2.n_bits:
        raise ValueError(f"Fingerprint bit lengths mismatch: {fp1.n_bits} vs {fp2.n_bits}")

    a = 0
    b = 0
    c = 0
    for b1, b2 in zip(fp1.bits, fp2.bits):
        if b1:
            a += 1
        if b2:
            b += 1
        if b1 and b2:
            c += 1
    return a, b, c


def tanimoto_similarity(fp1: MorganFingerprint, fp2: MorganFingerprint) -> float:
    """
    Computes Tanimoto (Jaccard) similarity coefficient:
    T(A, B) = c / (a + b - c)
    Range: [0.0, 1.0]
    """
    a, b, c = _get_intersection_counts(fp1, fp2)
    denom = a + b - c
    if denom == 0:
        return 1.0 if a == 0 and b == 0 else 0.0
    return c / denom


def dice_similarity(fp1: MorganFingerprint, fp2: MorganFingerprint) -> float:
    """
    Computes Sørensen-Dice similarity coefficient:
    D(A, B) = 2c / (a + b)
    Range: [0.0, 1.0]
    """
    a, b, c = _get_intersection_counts(fp1, fp2)
    denom = a + b
    if denom == 0:
        return 1.0 if a == 0 and b == 0 else 0.0
    return (2.0 * c) / denom


def cosine_similarity(fp1: MorganFingerprint, fp2: MorganFingerprint) -> float:
    """
    Computes Cosine (Otsuka-Ochiai) similarity coefficient:
    C(A, B) = c / sqrt(a * b)
    Range: [0.0, 1.0]
    """
    a, b, c = _get_intersection_counts(fp1, fp2)
    if a == 0 or b == 0:
        return 1.0 if a == 0 and b == 0 else 0.0
    return c / math.sqrt(a * b)


def tversky_index(
    fp1: MorganFingerprint,
    fp2: MorganFingerprint,
    alpha: float = 0.5,
    beta: float = 0.5,
) -> float:
    """
    Computes asymmetric Tversky similarity index:
    S(A, B; alpha, beta) = c / (c + alpha * (a - c) + beta * (b - c))
    Special cases:
    - alpha=1, beta=1 -> Tanimoto
    - alpha=0.5, beta=0.5 -> Dice
    - alpha=0.95, beta=0.05 -> Substructure search bias (prototypical similarity)
    """
    a, b, c = _get_intersection_counts(fp1, fp2)
    denom = c + alpha * (a - c) + beta * (b - c)
    if denom == 0:
        return 1.0 if a == 0 and b == 0 else 0.0
    return c / denom


def tanimoto_distance(fp1: MorganFingerprint, fp2: MorganFingerprint) -> float:
    """
    Computes Jaccard/Tanimoto metric distance:
    d_T(A, B) = 1.0 - T(A, B)
    Satisfies metric properties: non-negativity, symmetry, identity, and triangle inequality.
    """
    return 1.0 - tanimoto_similarity(fp1, fp2)


def compute_similarity_matrix(
    fps: List[MorganFingerprint],
    metric: str = "tanimoto",
) -> List[List[float]]:
    """
    Computes pairwise symmetric similarity matrix for a list of fingerprints.
    Supported metrics: 'tanimoto', 'dice', 'cosine'.
    """
    func_map = {
        "tanimoto": tanimoto_similarity,
        "dice": dice_similarity,
        "cosine": cosine_similarity,
    }
    sim_func = func_map.get(metric.lower())
    if not sim_func:
        raise ValueError(f"Unsupported metric: '{metric}'. Choose from {list(func_map.keys())}")

    n = len(fps)
    matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        matrix[i][i] = 1.0
        for j in range(i + 1, n):
            score = sim_func(fps[i], fps[j])
            matrix[i][j] = score
            matrix[j][i] = score

    return matrix
