import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as f:
	required = f.read().splitlines()

setuptools.setup(
	name="AmpliconComparison",
	version="0.0.1",
	author="",
	author_email="",
	description="Cycle comparison",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/pypa/sampleproject",
	project_urls={
		"Bug Tracker": "https://github.com/pypa/sampleproject/issues",
	},
	classifiers=[
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
	],
	install_requires=required,
	packages=setuptools.find_packages(),
	python_requires=">=3.7, <=3.11",
	entry_points={
        'console_scripts': [
            'AmpliconComparison = src:compare.py'
        	]
    	},
	zip_safe = False
)
