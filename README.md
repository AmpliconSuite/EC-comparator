# EC-comparator
Comparing cycle decompositions across technologies and methods.

## Table of Contents

- [Requirements](#requirements)
- [Getting started with pip](#getting-started-with-pip)
- [Usage](#usage)
- [Installation from source](#installation-from-source)
- [Output description](#output-description)
- [Key arguments](#key-arguments)
- [Help](#help)
- [Build and install (for developers)](#build-and-install-for-developers)
- [License](#license)
- [Related paper](#related-paper)
- [Contributors](#contributors)

## Requirements

EC-comparator requires:

- Python 3.10
- BedTools 2.31.1

BedTools is an external command-line dependency used by pybedtools. It is not installed by pip, so it must be installed separately.

Or, alternatively you can install these using conda/mamba ( assumes conda/mamba is preinstalled):

```
mamba install -c conda-forge -c bioconda python=3.10 bedtools=2.31.1
```


## Getting started with pip


``` 
python -m pip install EC-comparator==0.1.0
```

## Usage

Run the following example as test:

```bash
EC-comparator -a ./examples/ecdna1/true.bed \
               -b ./examples/ecdna1/reconstructed.bed \
               -d ./examples/ecdna1/output --plot --report
```

## Installation from source

Below is how to install `EC-comparator` from source. It assumes installed:
- conda/mamba
- git

```bash
git clone https://github.com/AmpliconSuite/EC-comparator.git
cd EC-comparator

mamba env create -f environment.yml
conda activate eccomparator

python -m pip install -r requirements.txt && python -m pip install .
which EC-comparator
EC-comparator --version
```

## Output description

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

| File | Description |
|---|---|
| `breakpoints_matched.txt` | Matched breakpoints between `s1` and `s2` |
| `breakpoints_profile_s1.txt` | Breakpoint profile of structure `s1` |
| `breakpoints_profile_s2.txt` | Breakpoint profile of structure `s2` |
| `coverage_breakpoints_profile.pdf` / `.png` / `.svg` | Coverage profile plot with breakpoints overlaid (same plot, three formats) |
| `coverage_profile.pdf` / `.png` / `.svg` | Coverage (copy-number) profile plot (same plot, three formats) |
| `coverage_profile_s1.txt` | Coverage profile values for structure `s1` |
| `coverage_profile_s2.txt` | Coverage profile values for structure `s2` |
| `metrics.json` | Computed distances and the configuration used to compute them |
| `report.html` / `report.pdf` | Generated comparison report (HTML and PDF versions) |
| `s1_input_filtered.bed` | Input structure `s1`, filtered (see `--min-cn`) |
| `s2_input_filtered.bed` | Input structure `s2`, filtered (see `--min-cn`) |
| `total_cost.png` | Radial plot of the distances / total cost |
| `total_cost_table.png` | Table image of the distances / total cost |

`coverage_breakpoints_profile.png`
![coverage_breakpoints_profile.png](./examples/ecdna1/output/coverage_breakpoints_profile.png)

`original_structures.png`
![original_structures.png](./examples/ecdna1/output/original_structures.png)

`total_cost`
| | |
|---|---|
| ![total_cost.png](./examples/ecdna1/output/total_cost.png) | ![total_cost_table.png](./examples/ecdna1/output/total_cost_table.png) |

`metrics.json`
```json
{
    "configs": {
        "breakpoint_dist": {
            ...
            },
            "gaussian": {
                "weight": 1,
                "threshold": 2,
                "enable": false,
                "sigma": 500,
                "amplitude": 1,
                "breakpoint_dist": {
                    "weight": 1,
                    "enable": true
                }
            },
            ...
        },
...
    },
    "distances": {
        "cn_hamming_dist": 350,
        "cn_hamming_norm_dist": 0.05,
        "cn_cos_dist": 0.01,
        "cn_jc_dist": 0.1,
        "breakpoint_dist": 0.64,
        "fragments_dist": 0.16,
        "cycles_dist": 0.05,
        "total_cost": 0.9,
       ...
    }
}%      
```

### Distances explained

Here we compute the cycle distance between `s1` and `s2`.

| ID | Distance | Description |
|---|---|---|
| d1 | Hamming distance | measures differences of genomic regions presence / absence for s1 and s2 |
| d2 | Cosine distance | measures the copy-number distance between the s1 and s2 using cosine distance |
| d3 | Min-max distance | measures the copy-number distance between the s1 and s2 using jaccard distance |
| d4 | Bins distance | measures differences of genomic bins usage |
| d5 | Paths distance | measures differences of paths decomposition between s1 and s2 |
| d6 | Breakpoint distance | measures differences between the breakpoint junctions profiles between s1 and s2 |

The radial plot shows all these distances, with values between 0 (low cost, high similarity) and 1 (high cost, low similarity). The more colorful the more distant are `s1` and `s2`.

**Total cost**: final score showing how dissimilar two reconstructions are (0 - highly similar, 5 - dissimilar)

## Key arguments

### Required arguments

| Flag | Value | Description |
|---|---|---|
| `-a`, `--first-structure` | `FIRST_STRUCTURE` | First structure (BED-like format) |
| `-b`, `--second-structure` | `SECOND_STRUCTURE` | Second structure (BED-like format) |
| `-d`, `--outdir` | `OUTDIR` | Output directory |

### Selected optional arguments

| Flag | Value | Description |
|---|---|---|
| `--plot` / `--no-plot` | — | Plot coverage profiles |
| `--report` / `--no-report` | — | Generate report (this the flag is set, it will also set 'plot') |
| `--gap` | `GAP` | Merge neighboring intervals within < gap (default: 1000000) |
| `--gaussian-sigma` | `GAUSSIAN_SIGMA` | Define standard deviation (default: 500) |

See [Help](#help) below for the full list of arguments.

## Help

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

## Build and install (for developers)

[Build and packaging](docs/dev.md)

## License

MIT

## Related paper

tba

## Contributors

- github@madagiurgiu25, Madalina Giurgiu-Kraljic, EMBL Heidelberg
- copilot, chatgpt, claude