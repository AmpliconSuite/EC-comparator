# EC-comparator
Comparing cycle decompositions across technologies and methods.

### Installation

Please install `python3.9` and run the following code to install all dependencies:

```bash
cd AmpliconComparison
mamba env create -f environment.yml
conda activate comparator
python -m pip install -r requirements.txt
```

And build and install (for developers):

```bash
python -m pip install build installer toml setuptools

cd EC-comparator
python -m build
python -m pip install --force-reinstall dist/EC-comparator-0.0.2-py3-none-any.whl
```

### Usage

Run the following test:

```bash
cd EC-comparator/AmpliconComparison
python main.py -a ../examples/ecdna1/true.bed \
               -b ../examples/ecdna1/reconstructed.bed \
               -d ../examples/ecdna1/output
```

### Output description

```bash
../examples/ecdna1/output/
|-- breakpoints_profile_s1.txt
|-- breakpoints_profile_s2.txt
|-- coverage_breakpoints_profile.png
|-- coverage_profile.png
|-- coverage_profile_s1.txt
|-- e1_coverage_profile.txt
|-- e2_coverage_profile.txt
|-- metrics.json
|-- report.html
|-- total_cost.png
|-- total_cost_table.png
```

### Help

```bash
usage: EC-comparator [-h] -a FIRST_STRUCTURE -b SECOND_STRUCTURE -d OUTDIR [--plot | --no-plot] [--report | --no-report] [--cn-hamming-dist CN_HAMMING_DIST] [--cn-cosine-dist CN_COSINE_DIST] [--cn-jc-dist CN_JC_DIST]
                          [--fragments-dist FRAGMENTS_DIST] [--cycles-dist CYCLES_DIST] [--breakpoint-dist BREAKPOINT_DIST]

EC-comparator - compare cycle sets

optional arguments:
  -h, --help            show this help message and exit
  -a FIRST_STRUCTURE, --first-structure FIRST_STRUCTURE
                        First structure (bed format)
  -b SECOND_STRUCTURE, --second-structure SECOND_STRUCTURE
                        Second structure (bed format)
  -d OUTDIR, --outdir OUTDIR
                        Output directory
  --plot, --no-plot     Plot coverage profiles
  --report, --no-report
                        Generate report (this will set 'plot' also on True)
  --cn-hamming-dist CN_HAMMING_DIST
                        Hamming distance between genomic footprint. Recommended when no copy-number information available (default: True)
  --cn-cosine-dist CN_COSINE_DIST
                        Cosine distance between the coverage profile.
  --cn-jc-dist CN_JC_DIST
                        Min-max distance / Jaccard distance between the coverage profile.
  --fragments-dist FRAGMENTS_DIST
                        Quantify the distance between fragments.
  --cycles-dist CYCLES_DIST
                        Quantify the distance between cycles.
  --breakpoint-dist BREAKPOINT_DIST
                        Quantify the distance between cycles.


```

### License

tbd

