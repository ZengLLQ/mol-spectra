"""
Molecular Graph representation with ring perception and aromaticity detection.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set, Tuple


class BondType(Enum):
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 4


# Standard IUPAC atomic weights (g/mol)
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
}

# Standard atomic numbers
ATOMIC_NUMBERS: Dict[str, int] = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Br": 35,
    "I": 53,
}

# Topological polar surface area contributions (Ertl et al., 2000)
TPSA_CONTRIBUTIONS: Dict[Tuple[str, int, int], float] = {
    # (symbol, valence, num_hydrogens) -> polar surface area (Angstrom^2)
    ("N", 3, 0): 3.24,
    ("N", 3, 1): 12.03,
    ("N", 3, 2): 26.02,
    ("N", 2, 0): 12.89,  # Pyridine-like
    ("N", 2, 1): 15.79,  # Pyrrole-like
    ("O", 2, 0): 9.23,   # Ether / ester
    ("O", 2, 1): 20.23,  # Hydroxyl
    ("O", 1, 0): 17.07,  # Carbonyl
    ("S", 2, 0): 25.30,  # Thioether
    ("S", 2, 1): 38.80,  # Thiol
    ("S", 1, 0): 32.50,  # Sulfoxide/sulfone O=S
    ("P", 3, 0): 0.0,
    ("P", 4, 0): 9.81,
}


@dataclass
class Atom:
    idx: int
    symbol: str
    formal_charge: int = 0
    implicit_h: int = 0
    is_aromatic: bool = False
    atomic_num: int = field(init=False)

    def __post_init__(self):
        self.atomic_num = ATOMIC_NUMBERS.get(self.symbol, 0)


@dataclass
class Bond:
    idx: int
    atom1_idx: int
    atom2_idx: int
    bond_type: BondType = BondType.SINGLE
    is_aromatic: bool = False

    def other_atom(self, atom_idx: int) -> int:
        if atom_idx == self.atom1_idx:
            return self.atom2_idx
        elif atom_idx == self.atom2_idx:
            return self.atom1_idx
        raise ValueError(f"Atom {atom_idx} is not part of bond {self.idx}")


class MolecularGraph:
    """
    Graph representation of a chemical molecule.
    Nodes are atoms; edges are covalent or aromatic bonds.
    """

    def __init__(self, name: str = "Molecule"):
        self.name = name
        self.atoms: Dict[int, Atom] = {}
        self.bonds: Dict[int, Bond] = {}
        self._adj: Dict[int, List[int]] = {}  # atom_idx -> list of bond_indices
        self._rings: Optional[List[List[int]]] = None

    def add_atom(
        self,
        symbol: str,
        formal_charge: int = 0,
        implicit_h: int = 0,
        is_aromatic: bool = False,
    ) -> int:
        idx = len(self.atoms)
        atom = Atom(
            idx=idx,
            symbol=symbol,
            formal_charge=formal_charge,
            implicit_h=implicit_h,
            is_aromatic=is_aromatic,
        )
        self.atoms[idx] = atom
        self._adj[idx] = []
        self._rings = None
        return idx

    def add_bond(
        self,
        atom1_idx: int,
        atom2_idx: int,
        bond_type: BondType = BondType.SINGLE,
        is_aromatic: bool = False,
    ) -> int:
        if atom1_idx not in self.atoms or atom2_idx not in self.atoms:
            raise KeyError("Both atoms must exist in the molecule before adding a bond.")
        if atom1_idx == atom2_idx:
            raise ValueError("Self-loop bonds are physically invalid.")

        bond_idx = len(self.bonds)
        bond = Bond(
            idx=bond_idx,
            atom1_idx=atom1_idx,
            atom2_idx=atom2_idx,
            bond_type=bond_type,
            is_aromatic=is_aromatic,
        )
        self.bonds[bond_idx] = bond
        self._adj[atom1_idx].append(bond_idx)
        self._adj[atom2_idx].append(bond_idx)
        self._rings = None
        return bond_idx

    def get_atom(self, atom_idx: int) -> Atom:
        return self.atoms[atom_idx]

    def get_bond(self, bond_idx: int) -> Bond:
        return self.bonds[bond_idx]

    def get_neighbors(self, atom_idx: int) -> List[int]:
        """Returns list of neighbor atom indices."""
        neighbors = []
        for b_idx in self._adj.get(atom_idx, []):
            bond = self.bonds[b_idx]
            neighbors.append(bond.other_atom(atom_idx))
        return neighbors

    def get_bond_between(self, atom1_idx: int, atom2_idx: int) -> Optional[Bond]:
        """Returns Bond object connecting atom1 and atom2 if exists."""
        for b_idx in self._adj.get(atom1_idx, []):
            bond = self.bonds[b_idx]
            if (bond.atom1_idx == atom1_idx and bond.atom2_idx == atom2_idx) or (
                bond.atom1_idx == atom2_idx and bond.atom2_idx == atom1_idx
            ):
                return bond
        return None

    def heavy_atom_count(self) -> int:
        """Counts non-hydrogen atoms."""
        return sum(1 for a in self.atoms.values() if a.symbol != "H")

    def total_atom_count(self) -> int:
        """Counts all atoms including implicit hydrogens."""
        explicit = len(self.atoms)
        implicit = sum(a.implicit_h for a in self.atoms.values())
        return explicit + implicit

    def molecular_weight(self) -> float:
        """Calculates molecular mass (g/mol) including implicit hydrogens."""
        mass = sum(ATOMIC_WEIGHTS.get(a.symbol, 0.0) for a in self.atoms.values())
        h_mass = ATOMIC_WEIGHTS["H"] * sum(a.implicit_h for a in self.atoms.values())
        return mass + h_mass

    def get_rings(self) -> List[List[int]]:
        """
        Detects Smallest Set of Smallest Rings (SSSR) using cycle basis DFS.
        Returns a list of cycles, each cycle being a list of atom indices.
        """
        if self._rings is not None:
            return self._rings

        cycles: List[List[int]] = []
        visited_nodes: Set[int] = set()
        parent_map: Dict[int, Optional[int]] = {}

        def dfs(u: int, p: Optional[int], path: List[int]):
            visited_nodes.add(u)
            parent_map[u] = p
            path.append(u)

            for v in self.get_neighbors(u):
                if v == p:
                    continue
                if v in visited_nodes:
                    # Cycle found
                    if v in path:
                        cycle_start_idx = path.index(v)
                        cycle = path[cycle_start_idx:]
                        # Normalize cycle to prevent duplicates
                        sorted_cycle = sorted(cycle)
                        if sorted_cycle not in [sorted(c) for c in cycles]:
                            cycles.append(cycle)
                else:
                    dfs(v, u, path)
            path.pop()

        for atom_idx in self.atoms:
            if atom_idx not in visited_nodes:
                dfs(atom_idx, None, [])

        self._rings = cycles
        return cycles

    def perceive_aromaticity(self):
        """
        Detects aromatic rings (5 and 6 membered rings with alternating double bonds
        or heteroatoms satisfying Hückel's 4n+2 pi-electron rule).
        Updates Atom.is_aromatic and Bond.is_aromatic.
        """
        rings = self.get_rings()
        for ring in rings:
            n_atoms = len(ring)
            if n_atoms not in (5, 6):
                continue

            # Check pi-electron contribution
            pi_electrons = 0
            is_candidate = True

            for atom_idx in ring:
                atom = self.atoms[atom_idx]
                # Check bonds in ring
                attached_bonds = [
                    self.get_bond_between(atom_idx, other)
                    for other in ring
                    if self.get_bond_between(atom_idx, other) is not None
                ]
                has_double_or_aromatic = any(
                    b.bond_type in (BondType.DOUBLE, BondType.AROMATIC) or b.is_aromatic
                    for b in attached_bonds
                )

                if atom.symbol == "C":
                    if has_double_or_aromatic or atom.is_aromatic or (len(attached_bonds) == 2 and atom.implicit_h == 1):
                        pi_electrons += 1
                    else:
                        is_candidate = False
                        break
                elif atom.symbol == "N":
                    if has_double_or_aromatic or atom.is_aromatic:
                        pi_electrons += 1  # Pyridine-like
                    else:
                        pi_electrons += 2  # Pyrrole-like (lone pair)
                elif atom.symbol == "O" or atom.symbol == "S":
                    pi_electrons += 2  # Furan/thiophene-like (lone pair)
                else:
                    is_candidate = False
                    break

            # Hückel's Rule: 4n + 2 pi-electrons (for 5/6 rings, 6 pi-electrons -> n=1)
            if is_candidate and pi_electrons == 6:
                for atom_idx in ring:
                    self.atoms[atom_idx].is_aromatic = True
                for i in range(len(ring)):
                    a1 = ring[i]
                    a2 = ring[(i + 1) % len(ring)]
                    bond = self.get_bond_between(a1, a2)
                    if bond:
                        bond.is_aromatic = True
                        bond.bond_type = BondType.AROMATIC

    def chemical_formula(self) -> str:
        """Computes Hill chemical formula string (C first, then H, then alphabetical)."""
        counts: Dict[str, int] = {}
        for a in self.atoms.values():
            counts[a.symbol] = counts.get(a.symbol, 0) + 1
            if a.implicit_h > 0:
                counts["H"] = counts.get("H", 0) + a.implicit_h

        formula_parts = []
        # C first
        if "C" in counts:
            c_cnt = counts.pop("C")
            formula_parts.append(f"C{c_cnt if c_cnt > 1 else ''}")
        # H second
        if "H" in counts:
            h_cnt = counts.pop("H")
            formula_parts.append(f"H{h_cnt if h_cnt > 1 else ''}")
        # Rest alphabetical
        for sym in sorted(counts.keys()):
            cnt = counts[sym]
            formula_parts.append(f"{sym}{cnt if cnt > 1 else ''}")

        return "".join(formula_parts)
