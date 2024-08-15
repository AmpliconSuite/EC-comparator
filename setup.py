import setuptools
from setuptools.command.install import install
import os
import shutil
import sys
import toml

class UninstallCommand(setuptools.Command):
	"""Custom command to uninstall the package."""

	user_options = []

	def initialize_options(self):
		pass

	def finalize_options(self):
		pass

	def run(self):
		# Uninstall the package
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
		# First, run the uninstall command
		self.run_command('uninstall')

		# Then proceed with the standard install
		install.run(self)

def parse_requirements(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Parse requirements
install_requires = parse_requirements('requirements.txt')

config = {}
# Load the configuration from pyproject.toml
with open("pyproject.toml", "r") as f:
	config = toml.load(f)

# with open("README.md", "r", encoding="utf-8") as fh:
# 	config["project"]["long_description"] = fh.read()
# 	# long_description = fh.read()


# Extract the 'project' section
project_config = config["project"]
print(project_config)

# Convert the authors to the format expected by setuptools
authors = project_config.get("authors", [])
author_names = [author.get("name") for author in authors]
author_emails = [author.get("email") for author in authors]

setuptools.setup(
	name=project_config["name"],
	version=project_config["version"],
	description=project_config.get("description", ""),
	long_description=open(project_config.get("readme", "README.md")).read(),
	long_description_content_type="text/markdown",
	author=", ".join(author_names),
	author_email=", ".join(author_emails),
	license=project_config.get("license", {}).get("text", ""),
	classifiers=project_config.get("classifiers", []),
	keywords=project_config.get("keywords", []),
	url=project_config.get("urls", {}).get("homepage", ""),
	packages=setuptools.find_packages(),
	install_requires=install_requires,
	# extras_require=project_config.get("optional-dependencies",{}).get("dev",{}),
	# name="AmpliconComparison",
	# version="0.0.2",
	# author="",
	# author_email="",
	# description="Cycle comparison",
	# long_description=long_description,
	# long_description_content_type="text/markdown",
	# url="https://github.com/pypa/sampleproject",
	project_urls={
		"Bug Tracker": "https://github.com/pypa/sampleproject/issues",
	},
	# classifiers=[
	# 	"License :: OSI Approved :: MIT License",
	# 	"Operating System :: OS Independent",
	# ],
	# scripts=["AmpliconComparison/main.py"],
	entry_points={
		'console_scripts': [
			'AmpliconComparison = AmpliconComparison.main:main'
			]
		},
	cmdclass={
		# 'install': CustomInstallCommand,
		'uninstall': UninstallCommand,
	},
	zip_safe = False
)
