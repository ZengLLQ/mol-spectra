"""
Bioavailability, QSAR, and Drug-Likeness Filters (Lipinski Rule of Five, Veber, TPSA, LogP).
"""

from dataclasses import dataclass, field
from typing import List, Tuple
from .graph import MolecularGraph, BondType, TPSA_CONTRIBUTIONS


@dataclass
class LipinskiProfile:
    molecular_weight: float
    log_p: float
    hbd: int  # Hydrogen Bond Donors
    hba: int  # Hydrogen Bond Acceptors
    violations: int
    passed: bool
    violation_reasons: List[str] = field(default_factory=list)


@dataclass
class VeberProfile:
    rotatable_bonds: int
    tpsa: float
    violations: int
    passed: bool
    violation_reasons: List[str] = field(default_factory=list)


def calculate_hbd(mol: MolecularGraph) -> int:
    """
    Counts Hydrogen Bond Donors (HBD).
    Number of OH and NH groups (nitrogen and oxygen atoms bonded to at least one hydrogen).
    """
    hbd = 0
    for atom in mol.atoms.values():
        if atom.symbol in ("N", "O"):
            # Check implicit hydrogens
            if atom.implicit_h > 0:
                hbd += atom.implicit_h
            else:
                # Check explicit H neighbors
                for n_idx in mol.get_neighbors(atom.idx):
                    if mol.atoms[n_idx].symbol == "H":
                        hbd += 1
    return hbd


def calculate_hba(mol: MolecularGraph) -> int:
    """
    Counts Hydrogen Bond Acceptors (HBA).
    Number of nitrogen and oxygen atoms.
    """
    return sum(1 for a in mol.atoms.values() if a.symbol in ("N", "O"))


def count_rotatable_bonds(mol: MolecularGraph) -> int:
    """
    Counts rotatable bonds according to Veber criteria:
    Any single non-ring bond, bounded to non-hydrogen atoms, where neither atom is terminal.
    Excludes amide C-N bonds due to partial double bond character.
    """
    mol.get_rings()
    ring_bonds = set()
    for ring in mol.get_rings():
        for i in range(len(ring)):
            a1 = ring[i]
            a2 = ring[(i + 1) % len(ring)]
            b = mol.get_bond_between(a1, a2)
            if b:
                ring_bonds.add(b.idx)

    rotatable = 0
    for bond in mol.bonds.values():
        if bond.idx in ring_bonds:
            continue
        if bond.bond_type != BondType.SINGLE:
            continue

        a1 = mol.atoms[bond.atom1_idx]
        a2 = mol.atoms[bond.atom2_idx]

        # Ignore bonds to hydrogen
        if a1.symbol == "H" or a2.symbol == "H":
            continue

        # Check heavy atom neighbors for both ends (non-terminal criteria)
        n1_heavy = sum(1 for n in mol.get_neighbors(a1.idx) if mol.atoms[n].symbol != "H")
        n2_heavy = sum(1 for n in mol.get_neighbors(a2.idx) if mol.atoms[n].symbol != "H")

        if n1_heavy > 1 and n2_heavy > 1:
            # Check for amide bond (C(=O)-N)
            is_amide = False
            for c_atom, n_atom in [(a1, a2), (a2, a1)]:
                if c_atom.symbol == "C" and n_atom.symbol == "N":
                    for cn in mol.get_neighbors(c_atom.idx):
                        if cn != n_atom.idx and mol.atoms[cn].symbol == "O":
                            cbond = mol.get_bond_between(c_atom.idx, cn)
                            if cbond and cbond.bond_type == BondType.DOUBLE:
                                is_amide = True
                                break
            if not is_amide:
                rotatable += 1

    return rotatable


def calculate_tpsa(mol: MolecularGraph) -> float:
    """
    Calculates Topological Polar Surface Area (TPSA in Angstrom^2) based on Ertl et al. (2000).
    Sums polar contributions from heteroatoms (N, O, S, P).
    """
    tpsa = 0.0
    for atom in mol.atoms.values():
        if atom.symbol not in ("N", "O", "S", "P"):
            continue

        valence = len(mol._adj.get(atom.idx, []))
        num_h = atom.implicit_h + sum(
            1 for n in mol.get_neighbors(atom.idx) if mol.atoms[n].symbol == "H"
        )

        # Lookup in contribution table or fallback approximation
        key = (atom.symbol, valence, num_h)
        contrib = TPSA_CONTRIBUTIONS.get(key)
        if contrib is None:
            # Generic fallback by atom and H count
            if atom.symbol == "O":
                contrib = 20.23 if num_h > 0 else 9.23
            elif atom.symbol == "N":
                contrib = 12.03 if num_h > 0 else 3.24
            elif atom.symbol == "S":
                contrib = 38.80 if num_h > 0 else 25.30
            else:
                contrib = 10.0
        tpsa += contrib

    return round(tpsa, 2)


def estimate_logp(mol: MolecularGraph) -> float:
    """
    Estimates octanol-water partition coefficient (logP) using atomic/fragment additivity
    (simplified Wildman-Crippen approach).
    """
    mol.perceive_aromaticity()
    logp = 0.0

    for atom in mol.atoms.values():
        sym = atom.symbol
        if sym == "C":
            if atom.is_aromatic:
                logp += 0.25
            else:
                logp += 0.10
        elif sym == "N":
            if atom.is_aromatic:
                logp -= 0.45
            else:
                logp -= 0.85
        elif sym == "O":
            if atom.implicit_h > 0:
                logp -= 0.40  # Hydroxyl
            else:
                logp -= 0.15  # Ether / ester
        elif sym == "F":
            logp += 0.38
        elif sym == "Cl":
            logp += 0.65
        elif sym == "Br":
            logp += 0.85
        elif sym == "I":
            logp += 1.10
        elif sym == "S":
            logp += 0.45
        elif sym == "P":
            logp += 0.25

        # Hydrogens attached to carbons increase lipophilicity
        if sym == "C" and atom.implicit_h > 0:
            logp += 0.10 * atom.implicit_h

    return round(logp, 2)


def evaluate_lipinski(mol: MolecularGraph) -> LipinskiProfile:
    """
    Evaluates Lipinski's Rule of Five for oral drug-likeness:
    1. Molecular Weight <= 500 Da
    2. LogP <= 5.0
    3. HBD <= 5
    4. HBA <= 10
    A molecule passes if it has <= 1 violation.
    """
    mw = mol.molecular_weight()
    logp = estimate_logp(mol)
    hbd = calculate_hbd(mol)
    hba = calculate_hba(mol)

    violations = []
    if mw > 500.0:
        violations.append(f"MW ({mw:.1f}) > 500 Da")
    if logp > 5.0:
        violations.append(f"LogP ({logp:.2f}) > 5.0")
    if hbd > 5:
        violations.append(f"HBD ({hbd}) > 5")
    if hba > 10:
        violations.append(f"HBA ({hba}) > 10")

    passed = len(violations) <= 1

    return LipinskiProfile(
        molecular_weight=round(mw, 2),
        log_p=logp,
        hbd=hbd,
        hba=hba,
        violations=len(violations),
        passed=passed,
        violation_reasons=violations,
    )


def evaluate_veber(mol: MolecularGraph) -> VeberProfile:
    """
    Evaluates Veber's oral bioavailability criteria:
    1. Rotatable bonds <= 10
    2. TPSA <= 140 Angstrom^2
    """
    rot = count_rotatable_bonds(mol)
    tpsa = calculate_tpsa(mol)

    violations = []
    if rot > 10:
        violations.append(f"Rotatable bonds ({rot}) > 10")
    if tpsa > 140.0:
        violations.append(f"TPSA ({tpsa:.1f}) > 140 A^2")

    passed = len(violations) == 0

    return VeberProfile(
        rotatable_bonds=rot,
        tpsa=tpsa,
        violations=len(violations),
        passed=passed,
        violation_reasons=violations,
    )
