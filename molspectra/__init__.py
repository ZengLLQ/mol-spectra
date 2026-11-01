"""
MolSpectra: Topological Molecular Fingerprinting, Chemical Similarity, and Bioavailability Profiling.
"""

from .graph import Atom, Bond, BondType, MolecularGraph
from .fingerprint import MorganFingerprint, generate_ecfp4, generate_ecfp6
from .similarity import (
    tanimoto_similarity,
    dice_similarity,
    cosine_similarity,
    tversky_index,
    compute_similarity_matrix,
)
from .bioavailability import (
    LipinskiProfile,
    VeberProfile,
    evaluate_lipinski,
    evaluate_veber,
    calculate_tpsa,
    estimate_logp,
)
from .library import get_reference_molecule, list_reference_molecules

__version__ = "0.1.0"
__all__ = [
    "Atom",
    "Bond",
    "BondType",
    "MolecularGraph",
    "MorganFingerprint",
    "generate_ecfp4",
    "generate_ecfp6",
    "tanimoto_similarity",
    "dice_similarity",
    "cosine_similarity",
    "tversky_index",
    "compute_similarity_matrix",
    "LipinskiProfile",
    "VeberProfile",
    "evaluate_lipinski",
    "evaluate_veber",
    "calculate_tpsa",
    "estimate_logp",
    "get_reference_molecule",
    "list_reference_molecules",
]
