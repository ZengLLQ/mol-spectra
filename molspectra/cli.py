"""
MolSpectra Command-Line Interface.
"""

import argparse
import sys
from typing import List

from .graph import MolecularGraph
from .fingerprint import generate_morgan_fingerprint, generate_ecfp4, generate_ecfp6
from .similarity import (
    tanimoto_similarity,
    dice_similarity,
    cosine_similarity,
    compute_similarity_matrix,
)
from .bioavailability import evaluate_lipinski, evaluate_veber
from .library import get_reference_molecule, list_reference_molecules


def cmd_list(args):
    print("Available Reference Molecules in MolSpectra Library:")
    for name in list_reference_molecules():
        mol = get_reference_molecule(name)
        print(f"  • {mol.name:<15} Formula: {mol.chemical_formula():<10} MW: {mol.molecular_weight():.2f} g/mol")


def cmd_info(args):
    mol = get_reference_molecule(args.molecule)
    lipinski = evaluate_lipinski(mol)
    veber = evaluate_veber(mol)

    print(f"=== Molecular Specification: {mol.name} ===")
    print(f"Chemical Formula     : {mol.chemical_formula()}")
    print(f"Molecular Weight     : {lipinski.molecular_weight:.2f} g/mol")
    print(f"Heavy Atom Count     : {mol.heavy_atom_count()}")
    print(f"Total Atom Count     : {mol.total_atom_count()}")
    print(f"Rings Detected       : {len(mol.get_rings())}")
    print("\n--- Lipinski's Rule of Five ---")
    print(f"LogP (lipophilicity) : {lipinski.log_p:.2f} (criteria: <= 5.0)")
    print(f"H-Bond Donors (HBD)  : {lipinski.hbd} (criteria: <= 5)")
    print(f"H-Bond Acceptors(HBA): {lipinski.hba} (criteria: <= 10)")
    print(f"Lipinski Violations  : {lipinski.violations} -> {'PASSED' if lipinski.passed else 'FAILED'}")
    if lipinski.violation_reasons:
        for r in lipinski.violation_reasons:
            print(f"  ! {r}")

    print("\n--- Veber Bioavailability Profile ---")
    print(f"Rotatable Bonds      : {veber.rotatable_bonds} (criteria: <= 10)")
    print(f"TPSA (polar surface) : {veber.tpsa:.2f} A^2 (criteria: <= 140 A^2)")
    print(f"Veber Violations     : {veber.violations} -> {'PASSED' if veber.passed else 'FAILED'}")


def cmd_fingerprint(args):
    mol = get_reference_molecule(args.molecule)
    fp = generate_morgan_fingerprint(mol, radius=args.radius, n_bits=args.bits)
    print(f"=== ECFP Fingerprint: {mol.name} ===")
    print(f"Radius               : {args.radius} (ECFP{args.radius * 2})")
    print(f"Bit Vector Length    : {args.bits}")
    print(f"Active Bits (On)     : {fp.popcount} / {args.bits} ({fp.popcount / args.bits * 100:.2f}%)")
    print(f"Unique Substructures : {len(fp.raw_features)}")
    print(f"Active Bit Indices   : {fp.on_bits[:20]}..." if len(fp.on_bits) > 20 else f"Active Bit Indices   : {fp.on_bits}")


def cmd_similarity(args):
    m1 = get_reference_molecule(args.mol1)
    m2 = get_reference_molecule(args.mol2)
    fp1 = generate_ecfp4(m1, n_bits=1024)
    fp2 = generate_ecfp4(m2, n_bits=1024)

    t_score = tanimoto_similarity(fp1, fp2)
    d_score = dice_similarity(fp1, fp2)
    c_score = cosine_similarity(fp1, fp2)

    print(f"=== Pairwise Chemical Similarity: {m1.name} vs {m2.name} ===")
    print(f"ECFP4 Tanimoto (Jaccard) : {t_score:.4f}")
    print(f"ECFP4 Sørensen-Dice      : {d_score:.4f}")
    print(f"ECFP4 Cosine (Otsuka)    : {c_score:.4f}")


def cmd_matrix(args):
    names = list_reference_molecules()
    mols = [get_reference_molecule(n) for n in names]
    fps = [generate_ecfp4(m, n_bits=1024) for m in mols]

    matrix = compute_similarity_matrix(fps, metric=args.metric)

    print(f"=== Cross-Similarity Matrix ({args.metric.capitalize()}, ECFP4) ===")
    header = f"{'':<14} | " + " | ".join(f"{n[:8]:<8}" for n in names)
    print(header)
    print("-" * len(header))
    for i, row in enumerate(matrix):
        row_str = f"{names[i]:<14} | " + " | ".join(f"{val:.4f}  " for val in row)
        print(row_str)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="molspectra",
        description="MolSpectra: Topological Molecular Fingerprints, Chemical Similarity, and Bioavailability",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # list
    subparsers.add_parser("list", help="List available reference molecules")

    # info
    p_info = subparsers.add_parser("info", help="Display molecular properties, Lipinski, and Veber profile")
    p_info.add_argument("molecule", type=str, help="Molecule name (e.g., aspirin, ibuprofen, paracetamol)")

    # fingerprint
    p_fp = subparsers.add_parser("fingerprint", help="Generate topological circular fingerprint")
    p_fp.add_argument("molecule", type=str, help="Molecule name")
    p_fp.add_argument("--radius", type=int, default=2, help="Neighborhood radius (default: 2 -> ECFP4)")
    p_fp.add_argument("--bits", type=int, default=1024, help="Bit vector length (default: 1024)")

    # similarity
    p_sim = subparsers.add_parser("similarity", help="Compute similarity between two molecules")
    p_sim.add_argument("mol1", type=str, help="First molecule")
    p_sim.add_argument("mol2", type=str, help="Second molecule")

    # matrix
    p_mat = subparsers.add_parser("matrix", help="Compute full pairwise similarity matrix")
    p_mat.add_argument("--metric", type=str, default="tanimoto", choices=["tanimoto", "dice", "cosine"])

    return parser


def main(argv: List[str] = None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "list":
        cmd_list(args)
    elif args.command == "info":
        cmd_info(args)
    elif args.command == "fingerprint":
        cmd_fingerprint(args)
    elif args.command == "similarity":
        cmd_similarity(args)
    elif args.command == "matrix":
        cmd_matrix(args)


if __name__ == "__main__":
    main()
