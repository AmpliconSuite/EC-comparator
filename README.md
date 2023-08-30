# AmpliconComparison
Comparing cycle decompositions across technologies and methods.

### Installation

Please install `python3` and run the following code to install all dependencies:

```bash
cd AmpliconComparison
python -m pip install .
```

### Usage

Run the following test:

```bash
cd AmpliconComparison/src
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
usage: AmpliconComparison [-h] -a FIRST_STRUCTURE -b SECOND_STRUCTURE -d OUTDIR [--cos-similarity COS_SIMILARITY] [--hamming-distance-norm HAMMING_DISTANCE_NORM]
                          [--cosine-distance COSINE_DISTANCE] [--fragments-overlap-norm FRAGMENTS_OVERLAP_NORM] [--cycles-overlap-norm CYCLES_OVERLAP_NORM]
                          [--euclidian-distance EUCLIDIAN_DISTANCE] [--euclidian-distance-threshold EUCLIDIAN_DISTANCE_THRESHOLD] [--relative-distance RELATIVE_DISTANCE]
                          [--relative-distance-threshold RELATIVE_DISTANCE_THRESHOLD]

AmpliconComparison - compare cycle sets

optional arguments:
  -h, --help            show this help message and exit
  -a FIRST_STRUCTURE, --first-structure FIRST_STRUCTURE
                        First structure (bed format)
  -b SECOND_STRUCTURE, --second-structure SECOND_STRUCTURE
                        Second structure (bed format)
  -d OUTDIR, --outdir OUTDIR
                        Output directory
  --cos-similarity COS_SIMILARITY
                        Cosine similarity between coverage tracks (default: True)
  --hamming-distance-norm HAMMING_DISTANCE_NORM
                        Hamming distance between genomic footprint. Recommended when no copy-number information available (default: True)
  --cosine-distance COSINE_DISTANCE
                        Cosine distance between the coverage profile.
  --fragments-overlap-norm FRAGMENTS_OVERLAP_NORM
                        Quantify the distance between fragments.
  --cycles-overlap-norm CYCLES_OVERLAP_NORM
                        Quantify the distance between cycles.
  --euclidian-distance EUCLIDIAN_DISTANCE
                        Use euclidian distance for breakpoint-pair matching (default: True)
  --euclidian-distance-threshold EUCLIDIAN_DISTANCE_THRESHOLD
                        Distance threshold to accept two breakpoint-pairs as matched (default: 1000)
  --relative-distance RELATIVE_DISTANCE
                        Relative distance score for breakpoint matching (default: False)
  --relative-distance-threshold RELATIVE_DISTANCE_THRESHOLD
                        Distance threshold to accept two breakpoint-pairs as matched (default: 0.3)

```

### License

tbd

