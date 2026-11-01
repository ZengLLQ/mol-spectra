"""
Morgan / Extended-Connectivity Fingerprint (ECFP) Generator.
Implements circular topological neighborhood hash aggregation and bit vector folding.
"""

import hashlib
import struct
from typing import Dict, List, Set, Tuple
from .graph import MolecularGraph, BondType


def _hash_int_sequence(seq: List[int]) -> int:
    """Hashes a sequence of integers into a deterministic 64-bit integer."""
    u64_seq = [x & 0xFFFFFFFFFFFFFFFF for x in seq]
    packed = struct.pack(f"{len(u64_seq)}Q", *u64_seq)
    digest = hashlib.sha256(packed).digest()
    # Take first 8 bytes as an unsigned 64-bit integer
    return struct.unpack("<Q", digest[:8])[0]


class MorganFingerprint:
    """
    Circular Extended-Connectivity Fingerprint representation.
    """

    def __init__(self, n_bits: int = 1024):
        self.n_bits = n_bits
        self.bits = [0] * n_bits
        self.raw_features: Set[int] = set()

    def set_bit(self, feature_hash: int):
        self.raw_features.add(feature_hash)
        bit_idx = feature_hash % self.n_bits
        self.bits[bit_idx] = 1

    @property
    def on_bits(self) -> List[int]:
        """Returns indices of all set (1) bits."""
        return [i for i, b in enumerate(self.bits) if b == 1]

    @property
    def popcount(self) -> int:
        """Returns total number of set bits (Hamming weight)."""
        return sum(self.bits)

    def to_bitstring(self) -> str:
        return "".join(str(b) for b in self.bits)

    def __repr__(self) -> str:
        return f"<MorganFingerprint bits={self.n_bits} on={self.popcount} features={len(self.raw_features)}>"


def generate_morgan_fingerprint(
    mol: MolecularGraph,
    radius: int = 2,
    n_bits: int = 1024,
) -> MorganFingerprint:
    """
    Generates circular topological Morgan / ECFP fingerprint.
    radius=2 corresponds to ECFP4 (diameter 4 bonds).
    radius=3 corresponds to ECFP6 (diameter 6 bonds).
    """
    # Ensure aromaticity is detected
    mol.perceive_aromaticity()

    fp = MorganFingerprint(n_bits=n_bits)
    if len(mol.atoms) == 0:
        return fp

    # Round 0: Calculate initial atom invariants
    current_invariants: Dict[int, int] = {}
    for atom_idx, atom in mol.atoms.items():
        neighbors = mol.get_neighbors(atom_idx)
        heavy_degree = sum(1 for n in neighbors if mol.atoms[n].symbol != "H")
        
        # Invariant components:
        # [atomic_number, heavy_degree, implicit_h, formal_charge, aromatic_flag]
        inv_vec = [
            atom.atomic_num,
            heavy_degree,
            atom.implicit_h,
            atom.formal_charge,
            1 if atom.is_aromatic else 0,
        ]
        h0 = _hash_int_sequence(inv_vec)
        current_invariants[atom_idx] = h0
        fp.set_bit(h0)

    # Iteration rounds 1 to radius
    for r in range(1, radius + 1):
        next_invariants: Dict[int, int] = {}

        for atom_idx in mol.atoms:
            # Collect neighbor tuples: (bond_type_code, neighbor_invariant)
            neighbor_features: List[Tuple[int, int]] = []
            for n_idx in mol.get_neighbors(atom_idx):
                bond = mol.get_bond_between(atom_idx, n_idx)
                b_code = 1
                if bond:
                    if bond.is_aromatic:
                        b_code = 4
                    elif bond.bond_type == BondType.DOUBLE:
                        b_code = 2
                    elif bond.bond_type == BondType.TRIPLE:
                        b_code = 3
                neighbor_features.append((b_code, current_invariants[n_idx]))

            # Lexicographically sort neighbor features for permutation invariance
            neighbor_features.sort()

            # Flatten into integer sequence for hashing:
            # [radius_round, center_invariant, b1, n1, b2, n2, ...]
            hasher_input = [r, current_invariants[atom_idx]]
            for b_code, n_inv in neighbor_features:
                hasher_input.append(b_code)
                hasher_input.append(n_inv)

            h_r = _hash_int_sequence(hasher_input)
            next_invariants[atom_idx] = h_r
            fp.set_bit(h_r)

        current_invariants = next_invariants

    return fp


def generate_ecfp4(mol: MolecularGraph, n_bits: int = 1024) -> MorganFingerprint:
    """Convenience alias for radius=2 ECFP4 fingerprint."""
    return generate_morgan_fingerprint(mol, radius=2, n_bits=n_bits)


def generate_ecfp6(mol: MolecularGraph, n_bits: int = 1024) -> MorganFingerprint:
    """Convenience alias for radius=3 ECFP6 fingerprint."""
    return generate_morgan_fingerprint(mol, radius=3, n_bits=n_bits)
