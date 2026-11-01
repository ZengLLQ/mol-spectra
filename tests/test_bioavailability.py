"""
Tests for Bioavailability Profiling (Lipinski, Veber, TPSA, LogP, Rotatable Bonds).
"""

from molspectra.library import (
    create_aspirin,
    create_ibuprofen,
    create_paracetamol,
    create_caffeine,
    create_benzene,
)
from molspectra.bioavailability import (
    calculate_hbd,
    calculate_hba,
    count_rotatable_bonds,
    calculate_tpsa,
    estimate_logp,
    evaluate_lipinski,
    evaluate_veber,
)


def test_aspirin_bioavailability():
    asp = create_aspirin()
    hbd = calculate_hbd(asp)
    hba = calculate_hba(asp)
    rot = count_rotatable_bonds(asp)
    tpsa = calculate_tpsa(asp)
    logp = estimate_logp(asp)

    assert hbd == 1  # -COOH OH
    assert hba == 4  # 4 oxygens
    assert rot == 3  # ring-C to COOH, ring-C to OAc, ester O-Ac
    assert tpsa > 50.0  # COOH + ester
    assert logp < 5.0

    lip = evaluate_lipinski(asp)
    assert lip.passed is True
    assert lip.violations == 0

    veb = evaluate_veber(asp)
    assert veb.passed is True
    assert veb.violations == 0


def test_ibuprofen_lipinski():
    ibu = create_ibuprofen()
    lip = evaluate_lipinski(ibu)
    assert lip.passed is True
    assert lip.violations == 0
    assert lip.hbd == 1
    assert lip.hba == 2


def test_paracetamol_lipinski():
    para = create_paracetamol()
    lip = evaluate_lipinski(para)
    assert lip.passed is True
    assert lip.violations == 0
    assert lip.hbd == 2  # Phenol -OH and amide -NH-
    assert lip.hba == 3  # N and 2 Os


def test_caffeine_profile():
    caff = create_caffeine()
    hbd = calculate_hbd(caff)
    hba = calculate_hba(caff)

    assert hbd == 0  # no OH or NH (all nitrogens are methyl-substituted or bridge)
    assert hba == 6  # 4 nitrogens + 2 oxygens

    lip = evaluate_lipinski(caff)
    assert lip.passed is True
    assert lip.violations == 0


def test_benzene_profile():
    bz = create_benzene()
    assert calculate_hbd(bz) == 0
    assert calculate_hba(bz) == 0
    assert calculate_tpsa(bz) == 0.0
    assert count_rotatable_bonds(bz) == 0
