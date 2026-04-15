# EC-comparator
Comparing cycle decompositions across technologies and methods.

### Installation

Please install `python3.9` and run the following code to install all dependencies:

```bash
cd AmpliconComparison
mamba env create -f environment.yml
conda activate comparator
python setup.py install
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
cd EC-comparator/eccomparator
python main.py -a ../examples/ecdna1/true.bed \
               -b ../examples/ecdna1/reconstructed.bed \
               -d ../examples/ecdna1/output
```

### Output description

```bash
../examples/ecdna1/output/
├── breakpoints_matched.txt
├── breakpoints_profile_s1.txt
├── breakpoints_profile_s2.txt
├── coverage_breakpoints_profile.pdf
├── coverage_breakpoints_profile.png
├── coverage_breakpoints_profile.svg
├── coverage_profile.pdf
├── coverage_profile.png
├── coverage_profile.svg
├── coverage_profile_s1.txt
├── coverage_profile_s2.txt
├── metrics.json
├── report.html
├── report.pdf
├── s1_input_filtered.bed
├── s2_input_filtered.bed
├── total_cost.png
├── total_cost_table.png
```

### Help

```bash
usage: EC-comparator [-h] -a FIRST_STRUCTURE -b SECOND_STRUCTURE -d OUTDIR [--plot | --no-plot] [--report | --no-report] [--min-cn MIN_CN] [--no-cn-hamming-dist] [--no-cn-cosine-dist] [--no-cn-jc-dist] [--no-fragments-dist] [--no-cycles-dist] [--no-breakpoint-dist]
                     [--breakpoint-dist-calc BREAKPOINT_DIST_CALC] [--gaussian-sigma GAUSSIAN_SIGMA] [--gaussian-amplitude GAUSSIAN_AMPLITUDE] [--breakpoint-cost-function BREAKPOINT_COST_FUNCTION] [--breakpoint-cost-function-threshold BREAKPOINT_COST_FUNCTION_THRESHOLD]
                     [--breakpoint-selection-distance BREAKPOINT_SELECTION_DISTANCE] [--breakpoint-selection-threshold BREAKPOINT_SELECTION_THRESHOLD] [--gap GAP] [--debug | --no-debug]

Method for comparing ecDNA structures (sets of cycle/paths).

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -a FIRST_STRUCTURE, --first-structure FIRST_STRUCTURE
                        First structure (BED-like format)
  -b SECOND_STRUCTURE, --second-structure SECOND_STRUCTURE
                        Second structure (BED-like format)
  -d OUTDIR, --outdir OUTDIR
                        Output directory

optional arguments:
  --plot, --no-plot     Plot coverage profiles
  --report, --no-report
                        Generate report (this the flag is set, it will also set 'plot')
  --min-cn MIN_CN       Minimal copy-number or coverage (filter out structures with a lower value then min-cn, default: 0)

optional arguments, which metrics to include:
  --no-cn-hamming-dist  Disable hamming distance between genomic footprint. Recommended when no copy-number information available.
  --no-cn-cosine-dist   Disable cosine distance between the coverage profile. Recommended when no copy-number information available.
  --no-cn-jc-dist       Disable min-max distance / Jaccard distance between the coverage profile. Recommended when no copy-number information available.
  --no-fragments-dist   Disable metric to quantify the distance between fragments.
  --no-cycles-dist      Disable metric to quantify the distance between cycles.
  --no-breakpoint-dist  Disable metric to quantify the distance between breakpoints.

optional arguments, fine tune breakpoint matching distance:
  --breakpoint-dist-calc BREAKPOINT_DIST_CALC
                        Define how to compute distance between breakpoints pairs (default: breakpoint_match_unweighted). Options: breakpoint_match_unweighted, breakpoint_match_cn_weighted, breakpoint_match_cn_weighted_avg, breakpoint_match_unweighted_confidence,
                        breakpoint_match_cn_weighted_confidence, breakpoint_match_cn_weighted_avg_confidence, breakpoint_gaussian_confidence_unweighted, breakpoint_gaussian_confidence_cn_weighted
  --gaussian-sigma GAUSSIAN_SIGMA
                        Define standard deviation (default: 500)
  --gaussian-amplitude GAUSSIAN_AMPLITUDE
                        Define amplitude of distribution (default: 1)
  --breakpoint-cost-function BREAKPOINT_COST_FUNCTION
                        Distance used to compute the cost matrix (default: euclidian). Options: euclidian, gaussian, match_score
  --breakpoint-cost-function-threshold BREAKPOINT_COST_FUNCTION_THRESHOLD
                        Distance used to compute the cost matrix (default: 3000)
  --breakpoint-selection-distance BREAKPOINT_SELECTION_DISTANCE
                        Distance for which two breakpoint-pairs are selected for matched candidates (default: manhattan). Options: manhattan, gaussian, match_score
  --breakpoint-selection-threshold BREAKPOINT_SELECTION_THRESHOLD
                        Threshold for which two breakpoint-pairs are considered matched candidates (default: 1000)
  --gap GAP             Merge neighboring intervals within < gap (default: 1000000)
  --debug, --no-debug   Debug structures (for developers)
```

### License

tbd

