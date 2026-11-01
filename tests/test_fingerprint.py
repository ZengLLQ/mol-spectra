"""
Tests for Morgan / ECFP Fingerprint generation and hashing properties.
"""

from molspectra.fingerprint import generate_morgan_fingerprint, generate_ecfp4, generate_ecfp6
from molspectra.library import create_aspirin, create_paracetamol, create_benzene


def test_fingerprint_generation():
    asp = create_aspirin()
    fp = generate_ecfp4(asp, n_bits=1024)

    assert fp.n_bits == 1024
    assert fp.popcount > 0
    assert len(fp.on_bits) == fp.popcount
    assert len(fp.raw_features) >= fp.popcount  # due to possible bit collisions


def test_fingerprint_radius_expansion():
    para = create_paracetamol()
    fp2 = generate_morgan_fingerprint(para, radius=1, n_bits=1024)
    fp4 = generate_ecfp4(para, n_bits=1024)
    fp6 = generate_ecfp6(para, n_bits=1024)

    # Increasing radius expands feature set
    assert len(fp2.raw_features) <= len(fp4.raw_features)
    assert len(fp4.raw_features) <= len(fp6.raw_features)


def test_bitstring_format():
    bz = create_benzene()
    fp = generate_ecfp4(bz, n_bits=512)
    bitstr = fp.to_bitstring()
    assert len(bitstr) == 512
    assert set(bitstr).issubset({"0", "1"})
    assert bitstr.count("1") == fp.popcount
