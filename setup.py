import re
from pathlib import Path
from setuptools import find_packages, setup

BASE_DIR = Path(__file__).resolve().parent

def get_version():
    init_file = BASE_DIR / "eccomparator" / "__init__.py"
    content = init_file.read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', content, re.MULTILINE)
    if not match:
        raise RuntimeError(f"Unable to find __version__ string in {init_file}")
    return match.group(1)

__version__ = get_version()

with open(BASE_DIR / "requirements.txt", encoding="utf-8") as f:
    requirements = [
        line.strip()
        for line in f
        if line.strip() and not line.startswith("#")
    ]

setup(
    name="EC-comparator",
    version=__version__,
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