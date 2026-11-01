"""
Tests for MolecularGraph, Atom, Bond, Ring perception, and Aromaticity.
"""

import pytest
from molspectra.graph import MolecularGraph, BondType, ATOMIC_WEIGHTS
from molspectra.library import create_benzene, create_aspirin, create_caffeine


def test_atom_and_bond_creation():
    mol = MolecularGraph("Water")
    o = mol.add_atom("O", implicit_h=2)
    assert mol.get_atom(o).symbol == "O"
    assert mol.get_atom(o).atomic_num == 8
    assert mol.heavy_atom_count() == 1
    assert mol.total_atom_count() == 3
    assert pytest.approx(mol.molecular_weight(), 0.01) == 18.015


def test_benzene_ring_and_aromaticity():
    benzene = create_benzene()
    rings = benzene.get_rings()
    assert len(rings) == 1
    assert len(rings[0]) == 6

    # Verify all atoms and bonds are perceived as aromatic
    for atom in benzene.atoms.values():
        assert atom.is_aromatic is True
    for bond in benzene.bonds.values():
        assert bond.is_aromatic is True
        assert bond.bond_type == BondType.AROMATIC

    assert benzene.chemical_formula() == "C6H6"
    assert pytest.approx(benzene.molecular_weight(), 0.1) == 78.11


def test_caffeine_fused_rings():
    caffeine = create_caffeine()
    rings = caffeine.get_rings()
    # Fused bicyclic system has 2 fundamental rings
    assert len(rings) >= 2
    assert caffeine.chemical_formula() == "C8H10N4O2"
    assert pytest.approx(caffeine.molecular_weight(), 0.2) == 194.19


def test_aspirin_structure():
    asp = create_aspirin()
    assert asp.heavy_atom_count() == 13
    assert asp.chemical_formula() == "C9H8O4"
    assert pytest.approx(asp.molecular_weight(), 0.2) == 180.16
