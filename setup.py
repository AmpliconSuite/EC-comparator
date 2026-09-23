from pathlib import Path

from setuptools import find_packages, setup


BASE_DIR = Path(__file__).resolve().parent

with open(BASE_DIR / "requirements.txt", encoding="utf-8") as f:
    requirements = [
        line.strip()
        for line in f
        if line.strip() and not line.startswith("#")
    ]

setup(
    name="EC-comparator",
    version="0.0.3",
    description="Compare ecDNA structures and amplicon sets across technologies and tools.",
    author="Madalina Giurgiu-Kraljic",
    url="https://github.com/AmpliconSuite/EC-comparator",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "EC-comparator=eccomparator.main:main",
        ],
    },
    python_requires=">=3.10",
)