"""
Tests for Chemical Similarity metrics (Tanimoto, Dice, Cosine, Tversky).
"""

import pytest
from molspectra.fingerprint import generate_ecfp4
from molspectra.similarity import (
    tanimoto_similarity,
    dice_similarity,
    cosine_similarity,
    tversky_index,
    tanimoto_distance,
    compute_similarity_matrix,
)
from molspectra.library import create_aspirin, create_ibuprofen, create_paracetamol


def test_tanimoto_identity_and_symmetry():
    asp = create_aspirin()
    ibu = create_ibuprofen()
    fp_a = generate_ecfp4(asp, n_bits=1024)
    fp_b = generate_ecfp4(ibu, n_bits=1024)

    # Identity: T(A, A) == 1.0
    assert pytest.approx(tanimoto_similarity(fp_a, fp_a), 1e-6) == 1.0
    assert pytest.approx(dice_similarity(fp_a, fp_a), 1e-6) == 1.0
    assert pytest.approx(cosine_similarity(fp_a, fp_a), 1e-6) == 1.0

    # Symmetry: T(A, B) == T(B, A)
    t_ab = tanimoto_similarity(fp_a, fp_b)
    t_ba = tanimoto_similarity(fp_b, fp_a)
    assert pytest.approx(t_ab, 1e-6) == t_ba

    # Distance identity and symmetry
    assert pytest.approx(tanimoto_distance(fp_a, fp_a), 1e-6) == 0.0
    assert pytest.approx(tanimoto_distance(fp_a, fp_b), 1e-6) == 1.0 - t_ab


def test_metric_triangle_inequality():
    asp = create_aspirin()
    ibu = create_ibuprofen()
    para = create_paracetamol()

    fp1 = generate_ecfp4(asp, n_bits=1024)
    fp2 = generate_ecfp4(ibu, n_bits=1024)
    fp3 = generate_ecfp4(para, n_bits=1024)

    d12 = tanimoto_distance(fp1, fp2)
    d23 = tanimoto_distance(fp2, fp3)
    d13 = tanimoto_distance(fp1, fp3)

    # Triangle inequality: d(A, C) <= d(A, B) + d(B, C)
    assert d13 <= d12 + d23 + 1e-9


def test_tversky_properties():
    asp = create_aspirin()
    ibu = create_ibuprofen()
    fp_a = generate_ecfp4(asp, n_bits=1024)
    fp_b = generate_ecfp4(ibu, n_bits=1024)

    # Tversky(1, 1) == Tanimoto
    tv_11 = tversky_index(fp_a, fp_b, alpha=1.0, beta=1.0)
    assert pytest.approx(tv_11, 1e-6) == tanimoto_similarity(fp_a, fp_b)

    # Tversky(0.5, 0.5) == Dice
    tv_half = tversky_index(fp_a, fp_b, alpha=0.5, beta=0.5)
    assert pytest.approx(tv_half, 1e-6) == dice_similarity(fp_a, fp_b)


def test_similarity_matrix():
    asp = create_aspirin()
    para = create_paracetamol()
    fps = [generate_ecfp4(asp), generate_ecfp4(para)]

    mat = compute_similarity_matrix(fps, metric="tanimoto")
    assert len(mat) == 2
    assert len(mat[0]) == 2
    assert pytest.approx(mat[0][0], 1e-6) == 1.0
    assert pytest.approx(mat[1][1], 1e-6) == 1.0
    assert pytest.approx(mat[0][1], 1e-6) == mat[1][0]
