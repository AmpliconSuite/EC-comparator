import setuptools
from setuptools.command.install import install
import os
import shutil
import sys

class UninstallCommand(setuptools.Command):
    """Custom command to uninstall the package."""
    user_options = []
    def initialize_options(self): pass
    def finalize_options(self): pass
    def run(self):
        dist_name = self.distribution.get_name()
        installed_path = next(p for p in sys.path if os.path.exists(os.path.join(p, dist_name + '.egg-info')))
        if installed_path:
            confirm = input(f"Do you want to uninstall the existing version of {dist_name}? [y/N]: ").lower()
            if confirm == 'y':
                egg_info_dir = os.path.join(installed_path, dist_name + '.egg-info')
                if os.path.isdir(egg_info_dir):
                    shutil.rmtree(egg_info_dir)
                    print(f"Removed {egg_info_dir}")
                package_dir = os.path.join(installed_path, dist_name)
                if os.path.isdir(package_dir):
                    shutil.rmtree(package_dir)
                    print(f"Removed {package_dir}")
            else:
                print("Uninstallation aborted.")
        else:
            print(f"{dist_name} is not installed.")

class CustomInstallCommand(install):
    """Custom install command that first uninstalls the package."""
    def run(self):
        self.run_command('uninstall')
        install.run(self)

def parse_requirements(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

install_requires = parse_requirements('requirements.txt')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="EC-comparator",
    version="0.0.2",
    description="Comparing cycle decompositions across technologies and methods.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="",
    author_email="test@example.com",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    keywords=["example", "package"],
    url="https://github.com/AmpliconSuite/EC-comparator",
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    python_requires=">=3.9",
    project_urls={
        "Bug Tracker": "https://github.com/AmpliconSuite/EC-comparator/issues",
        "Source Code": "https://github.com/AmpliconSuite/EC-comparator",
    },
    entry_points={
        'console_scripts': [
            'EC-comparator = eccomparator.main:main'
        ]
    },
    cmdclass={
        'uninstall': UninstallCommand,
    },
    zip_safe=False
)