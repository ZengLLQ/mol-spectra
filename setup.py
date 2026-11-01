from setuptools import setup, find_packages

setup(
    name="molspectra",
    version="0.1.0",
    description="Topological Molecular Fingerprints, Chemical Similarity, and Bioavailability Profiling",
    author="ZengLLQ",
    packages=find_packages(),
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "molspectra = molspectra.cli:main",
        ],
    },
)
