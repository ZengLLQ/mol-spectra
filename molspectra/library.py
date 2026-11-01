"""
Curated library of bioactive pharmaceutical molecules and chemical benchmarks.
"""

from typing import Dict, List
from .graph import MolecularGraph, BondType


def create_benzene() -> MolecularGraph:
    """Benzene (C6H6) - 6-membered aromatic ring."""
    m = MolecularGraph("Benzene")
    for _ in range(6):
        m.add_atom("C", implicit_h=1)
    for i in range(6):
        m.add_bond(i, (i + 1) % 6, BondType.AROMATIC, is_aromatic=True)
    m.perceive_aromaticity()
    return m


def create_aspirin() -> MolecularGraph:
    """Aspirin / Acetylsalicylic acid (C9H8O4)."""
    m = MolecularGraph("Aspirin")
    # Benzene ring carbons (0..5)
    c0 = m.add_atom("C", implicit_h=0)  # C bonded to COOH
    c1 = m.add_atom("C", implicit_h=0)  # C bonded to O-acetyl
    c2 = m.add_atom("C", implicit_h=1)
    c3 = m.add_atom("C", implicit_h=1)
    c4 = m.add_atom("C", implicit_h=1)
    c5 = m.add_atom("C", implicit_h=1)

    # Benzene ring bonds
    m.add_bond(c0, c1, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c1, c2, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c2, c3, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c3, c4, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c4, c5, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c5, c0, BondType.AROMATIC, is_aromatic=True)

    # Carboxylic acid (-COOH) at c0
    c_carb = m.add_atom("C", implicit_h=0)
    o_double1 = m.add_atom("O", implicit_h=0)
    o_single1 = m.add_atom("O", implicit_h=1)  # -OH donor
    m.add_bond(c0, c_carb, BondType.SINGLE)
    m.add_bond(c_carb, o_double1, BondType.DOUBLE)
    m.add_bond(c_carb, o_single1, BondType.SINGLE)

    # Acetoxy group (-OCOCH3) at c1
    o_ester = m.add_atom("O", implicit_h=0)
    c_ac = m.add_atom("C", implicit_h=0)
    o_double2 = m.add_atom("O", implicit_h=0)
    c_me = m.add_atom("C", implicit_h=3)
    m.add_bond(c1, o_ester, BondType.SINGLE)
    m.add_bond(o_ester, c_ac, BondType.SINGLE)
    m.add_bond(c_ac, o_double2, BondType.DOUBLE)
    m.add_bond(c_ac, c_me, BondType.SINGLE)

    m.perceive_aromaticity()
    return m


def create_paracetamol() -> MolecularGraph:
    """Paracetamol / Acetaminophen (C8H9NO2)."""
    m = MolecularGraph("Paracetamol")
    # Benzene ring (0..5)
    c0 = m.add_atom("C", implicit_h=0)  # C-OH
    c1 = m.add_atom("C", implicit_h=1)
    c2 = m.add_atom("C", implicit_h=1)
    c3 = m.add_atom("C", implicit_h=0)  # C-NH-Ac
    c4 = m.add_atom("C", implicit_h=1)
    c5 = m.add_atom("C", implicit_h=1)

    m.add_bond(c0, c1, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c1, c2, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c2, c3, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c3, c4, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c4, c5, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c5, c0, BondType.AROMATIC, is_aromatic=True)

    # Hydroxyl (-OH) at c0
    o_oh = m.add_atom("O", implicit_h=1)
    m.add_bond(c0, o_oh, BondType.SINGLE)

    # Acetamide (-NH-CO-CH3) at c3
    n_nh = m.add_atom("N", implicit_h=1)
    c_ac = m.add_atom("C", implicit_h=0)
    o_carbonyl = m.add_atom("O", implicit_h=0)
    c_me = m.add_atom("C", implicit_h=3)
    m.add_bond(c3, n_nh, BondType.SINGLE)
    m.add_bond(n_nh, c_ac, BondType.SINGLE)
    m.add_bond(c_ac, o_carbonyl, BondType.DOUBLE)
    m.add_bond(c_ac, c_me, BondType.SINGLE)

    m.perceive_aromaticity()
    return m


def create_ibuprofen() -> MolecularGraph:
    """Ibuprofen (C13H18O2)."""
    m = MolecularGraph("Ibuprofen")
    # Benzene ring (0..5)
    c0 = m.add_atom("C", implicit_h=0)  # C-isobutyl
    c1 = m.add_atom("C", implicit_h=1)
    c2 = m.add_atom("C", implicit_h=1)
    c3 = m.add_atom("C", implicit_h=0)  # C-CH(CH3)-COOH
    c4 = m.add_atom("C", implicit_h=1)
    c5 = m.add_atom("C", implicit_h=1)

    m.add_bond(c0, c1, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c1, c2, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c2, c3, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c3, c4, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c4, c5, BondType.AROMATIC, is_aromatic=True)
    m.add_bond(c5, c0, BondType.AROMATIC, is_aromatic=True)

    # Isobutyl at c0: -CH2-CH(CH3)2
    c_ch2 = m.add_atom("C", implicit_h=2)
    c_ch = m.add_atom("C", implicit_h=1)
    c_me1 = m.add_atom("C", implicit_h=3)
    c_me2 = m.add_atom("C", implicit_h=3)
    m.add_bond(c0, c_ch2, BondType.SINGLE)
    m.add_bond(c_ch2, c_ch, BondType.SINGLE)
    m.add_bond(c_ch, c_me1, BondType.SINGLE)
    m.add_bond(c_ch, c_me2, BondType.SINGLE)

    # Propionic acid at c3: -CH(CH3)-COOH
    c_alpha = m.add_atom("C", implicit_h=1)
    c_me3 = m.add_atom("C", implicit_h=3)
    c_acid = m.add_atom("C", implicit_h=0)
    o_double = m.add_atom("O", implicit_h=0)
    o_single = m.add_atom("O", implicit_h=1)

    m.add_bond(c3, c_alpha, BondType.SINGLE)
    m.add_bond(c_alpha, c_me3, BondType.SINGLE)
    m.add_bond(c_alpha, c_acid, BondType.SINGLE)
    m.add_bond(c_acid, o_double, BondType.DOUBLE)
    m.add_bond(c_acid, o_single, BondType.SINGLE)

    m.perceive_aromaticity()
    return m


def create_caffeine() -> MolecularGraph:
    """Caffeine (C8H10N4O2) - Purine xanthine alkaloid."""
    m = MolecularGraph("Caffeine")
    # 6-membered ring atoms (0..5): C1, C2, N3, C4, C5, N6
    c0 = m.add_atom("C", implicit_h=0)  # C=O
    o0 = m.add_atom("O", implicit_h=0)
    n1 = m.add_atom("N", implicit_h=0)  # N-CH3
    c_me1 = m.add_atom("C", implicit_h=3)
    c2 = m.add_atom("C", implicit_h=0)  # C=O
    o2 = m.add_atom("O", implicit_h=0)
    n3 = m.add_atom("N", implicit_h=0)  # N-CH3
    c_me2 = m.add_atom("C", implicit_h=3)
    c4 = m.add_atom("C", implicit_h=0)  # bridge
    c5 = m.add_atom("C", implicit_h=0)  # bridge

    # 6-ring bonds
    m.add_bond(c0, o0, BondType.DOUBLE)
    m.add_bond(c0, n1, BondType.SINGLE)
    m.add_bond(n1, c_me1, BondType.SINGLE)
    m.add_bond(n1, c2, BondType.SINGLE)
    m.add_bond(c2, o2, BondType.DOUBLE)
    m.add_bond(c2, n3, BondType.SINGLE)
    m.add_bond(n3, c_me2, BondType.SINGLE)
    m.add_bond(n3, c4, BondType.SINGLE)
    m.add_bond(c4, c5, BondType.DOUBLE)
    m.add_bond(c5, c0, BondType.SINGLE)

    # Imidazole 5-ring attached to c4, c5
    n7 = m.add_atom("N", implicit_h=0)  # N-CH3
    c_me3 = m.add_atom("C", implicit_h=3)
    c8 = m.add_atom("C", implicit_h=1)  # C-H
    n9 = m.add_atom("N", implicit_h=0)

    m.add_bond(c5, n7, BondType.SINGLE)
    m.add_bond(n7, c_me3, BondType.SINGLE)
    m.add_bond(n7, c8, BondType.SINGLE)
    m.add_bond(c8, n9, BondType.DOUBLE)
    m.add_bond(n9, c4, BondType.SINGLE)

    m.perceive_aromaticity()
    return m


_LIBRARY: Dict[str, MolecularGraph] = {
    "benzene": create_benzene,
    "aspirin": create_aspirin,
    "paracetamol": create_paracetamol,
    "ibuprofen": create_ibuprofen,
    "caffeine": create_caffeine,
}


def get_reference_molecule(name: str) -> MolecularGraph:
    """Retrieves reference molecule graph by name."""
    key = name.strip().lower()
    if key not in _LIBRARY:
        raise KeyError(f"Unknown reference molecule: '{name}'. Available: {list(_LIBRARY.keys())}")
    return _LIBRARY[key]()


def list_reference_molecules() -> List[str]:
    """Returns names of all available reference molecules."""
    return list(_LIBRARY.keys())
