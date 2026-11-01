"""
Tests for MolSpectra CLI commands.
"""

from molspectra.cli import main


def test_cli_list(capsys):
    main(["list"])
    captured = capsys.readouterr()
    assert "Aspirin" in captured.out
    assert "Ibuprofen" in captured.out
    assert "Caffeine" in captured.out


def test_cli_info(capsys):
    main(["info", "aspirin"])
    captured = capsys.readouterr()
    assert "Molecular Specification: Aspirin" in captured.out
    assert "Formula: C9H8O4" in captured.out or "C9H8O4" in captured.out
    assert "Rule of Five" in captured.out
    assert "PASSED" in captured.out


def test_cli_fingerprint(capsys):
    main(["fingerprint", "aspirin", "--radius", "2", "--bits", "512"])
    captured = capsys.readouterr()
    assert "ECFP Fingerprint: Aspirin" in captured.out
    assert "Bit Vector Length    : 512" in captured.out
    assert "Active Bits (On)" in captured.out


def test_cli_similarity(capsys):
    main(["similarity", "aspirin", "paracetamol"])
    captured = capsys.readouterr()
    assert "Tanimoto" in captured.out
    assert "Dice" in captured.out
    assert "Cosine" in captured.out


def test_cli_matrix(capsys):
    main(["matrix", "--metric", "tanimoto"])
    captured = capsys.readouterr()
    assert "Cross-Similarity Matrix" in captured.out
    assert "aspirin" in captured.out
    assert "1.0000" in captured.out
