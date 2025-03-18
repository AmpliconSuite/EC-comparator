"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 10:30 AM 07/28/23

AmpliconComparison entry point
"""

import argparse
import os
import traceback
import datetime
import pprint
import numpy as np

from AmpliconComparison.configs import Configs
from AmpliconComparison.utils.utils import DDT as d
from AmpliconComparison.utils.utils import HEADER as h
from AmpliconComparison.metrics import compare

# from AmpliconComparison.utils.utils import DDT as d
# from AmpliconComparison.utils.utils import HEADER as h


def format_print(dict):
    """
    Print configuration
    """
    print()
    # pprint.pprint(dict)
    print(
        f"### Candidates for breakpoint maching: {dict[h.BREAKPOINT_DISTANCE][d.UNMATCHED_DISTANCE]}, {dict[h.BREAKPOINT_DISTANCE][d.UNMATCHED_THRESHOLD]}"
    )
    print(
        f"### Cost matrix function for minimum matching bipartite graph: {dict[h.BREAKPOINT_DISTANCE][h.DEFAULT]}"
    )
    print(
        f"### Compute distance similarity using: {dict[h.BREAKPOINT_DISTANCE][h.BREAKPOINT_DISTANCE_CALCULATION]}"
    )
    print()


def config(args):
    """
    Configure comparison based on user spefication
    """
    dict = Configs.DICT_DIST

    # set all defaults
    d.GAUSSIAN_SIGMA = args.gaussian_sigma
    d.GAUSSIAN_AMPL = args.gaussian_amplitude
    d.BP_MATCH_UNWEIGHTED = args.breakpoint_dist_calc

    # genomic footprint distance settings
    dict[h.GENOMIC_FOOTPRINT][d.HAMMING_NORM][h.ENABLE] = (
        False if args.no_cn_hamming_dist else True
    )
    dict[h.GENOMIC_FOOTPRINT][d.COSINE_DISTANCE][h.ENABLE] = (
        False if args.no_cn_cosine_dist else True
    )
    dict[h.GENOMIC_FOOTPRINT][d.COPYNUMBER_JC][h.ENABLE] = (
        False if args.no_cn_jc_dist else True
    )
    dict[h.GENOMIC_FOOTPRINT][d.FRAGMENTS_DISTANCE][h.ENABLE] = (
        False if args.no_fragments_dist else True
    )
    dict[h.GENOMIC_FOOTPRINT][d.CYCLES_DISTANCE][h.ENABLE] = (
        False if args.no_cycles_dist else True
    )

    # how to compute distance between breakpoints
    dict[h.BREAKPOINT_DISTANCE][
        h.BREAKPOINT_DISTANCE_CALCULATION
    ] = args.breakpoint_dist_calc
    # if the breakpoint distance calculation includes confidence then the breakpoint selection methods needs to be gaussian
    if dict[h.BREAKPOINT_DISTANCE][h.BREAKPOINT_DISTANCE_CALCULATION] in [
        d.BP_MATCH_UNWEIGHTED_CONFIDENCE,
        d.BP_MATCH_CN_WEIGHTED_CONFIDENCE,
        d.BP_MATCH_CN_WEIGHTED_AVG_CONFIDENCE,
        d.BP_GAUSSIAN_CONFIDENCE,
    ]:
        args.breakpoint_selection_distance = d.GAUSSIAN
    # dict[h.BREAKPOINT_DISTANCE][d.MATCH_NONLINEAR] = args.match_breakpoints_nonlinear

    dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = False
    dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.THRESHOLD] = 3000

    dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = False
    dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.THRESHOLD] = 0.3

    dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][d.GAUSSIAN_DEF_SIGMA] = args.gaussian_sigma
    dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][
        d.GAUSSIAN_DEF_AMPL
    ] = args.gaussian_amplitude
    dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][h.ENABLE] = False
    dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][h.THRESHOLD] = 2 * args.gaussian_amplitude

    dict[h.BREAKPOINT_DISTANCE][
        d.UNMATCHED_DISTANCE
    ] = args.breakpoint_selection_distance
    if args.breakpoint_selection_distance == d.GAUSSIAN:
        dict[h.BREAKPOINT_DISTANCE][d.UNMATCHED_THRESHOLD] = 2 * args.gaussian_amplitude
    else:
        dict[h.BREAKPOINT_DISTANCE][
            d.UNMATCHED_THRESHOLD
        ] = args.breakpoint_selection_threshold

    # what to use for computing the cost function
    if args.breakpoint_cost_function == d.EUCLIDIAN:
        dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.EUCLIDIAN
        dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = True
        dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][
            h.THRESHOLD
        ] = args.breakpoint_cost_function_threshold
    elif args.breakpoint_cost_function == d.RELATIVE_METRIC:
        dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.RELATIVE_METRIC
        dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = True
        dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][
            h.THRESHOLD
        ] = args.breakpoint_cost_function_threshold
    elif args.breakpoint_cost_function == d.GAUSSIAN:
        dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.GAUSSIAN
        dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][h.ENABLE] = True
        dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][h.THRESHOLD] = (
            2 * args.gaussian_amplitude
        )

    # if breakpoint distance is disabled
    if args.no_breakpoint_dist:
        dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = None
        dict[h.BREAKPOINT_DISTANCE][d.GAUSSIAN][h.ENABLE] = False
        dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = False
        dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = False
        dict[h.BREAKPOINT_DISTANCE][h.BREAKPOINT_DISTANCE_CALCULATION] = None

    format_print(dict)
    return dict


def main():
    try:
        parser = argparse.ArgumentParser(
            prog="EC-comparator",
            description="Method for comparing ecDNA structures (sets of cycle/paths).",
            add_help=True,
        )
        # required
        required_args = parser.add_argument_group("required arguments")
        required_args.add_argument(
            "-a",
            "--first-structure",
            help="First structure (BED-like format)",
            required=True,
        )
        required_args.add_argument(
            "-b",
            "--second-structure",
            help="Second structure (BED-like format)",
            required=True,
        )
        required_args.add_argument(
            "-d", "--outdir", help="Output directory", required=True
        )

        # optional - plotting and filtering threshold
        optional_args = parser.add_argument_group("optional arguments")
        optional_args.add_argument(
            "--plot",
            help="Plot coverage profiles",
            action=argparse.BooleanOptionalAction,
        )
        optional_args.add_argument(
            "--report",
            help="Generate report (this the flag is set, it will also set 'plot')",
            action=argparse.BooleanOptionalAction,
        )
        optional_args.add_argument(
            "--min-cn",
            help="Minimal copy-number or coverage (filter out structures with a lower value then min-cn, default: %(default)s)",
            required=False,
            default=0,
            type=float,
        )

        # optional - which metrics to include
        optional_args_dist = parser.add_argument_group(
            "optional arguments, which metrics to include"
        )
        optional_args_dist.add_argument(
            "--no-cn-hamming-dist",
            help="Disable hamming distance between genomic footprint. Recommended when no copy-number information available.",
            required=False,
            action="store_true",
        )
        optional_args_dist.add_argument(
            "--no-cn-cosine-dist",
            help="Disable cosine distance between the coverage profile. Recommended when no copy-number information available.",
            required=False,
            action="store_true",
        )
        optional_args_dist.add_argument(
            "--no-cn-jc-dist",
            help="Disable min-max distance / Jaccard distance between the coverage profile. Recommended when no copy-number information available.",
            required=False,
            action="store_true",
        )
        optional_args_dist.add_argument(
            "--no-fragments-dist",
            help="Disable metric to quantify the distance between fragments. ",
            required=False,
            action="store_true",
        )
        optional_args_dist.add_argument(
            "--no-cycles-dist",
            help="Disable metric to quantify the distance between cycles.",
            required=False,
            action="store_true",
        )
        optional_args_dist.add_argument(
            "--no-breakpoint-dist",
            help="Disable metric to quantify the distance between breakpoints.",
            required=False,
            action="store_true",
        )

        # optional - configure breakpoint matching distance
        optional_args_bp = parser.add_argument_group(
            "optional arguments, fine tune breakpoint matching distance"
        )
        optional_args_bp.add_argument(
            "--breakpoint-dist-calc",
            help="Define how to compute distance between breakpoints pairs (default: %(default)s). Options: "
            + ", ".join(d.BP_OPTIONS),
            required=False,
            default=d.BP_MATCH_UNWEIGHTED,
            type=str,
        )

        optional_args_bp.add_argument(
            "--gaussian-sigma",
            help="Define standard deviation (default: %(default)s)",
            required=False,
            default=d.GAUSSIAN_SIGMA,
            type=int,
        )
        optional_args_bp.add_argument(
            "--gaussian-amplitude",
            help="Define amplitude of distribution (default: %(default)s)",
            required=False,
            default=d.GAUSSIAN_AMPL,
            type=int,
        )

        optional_args_bp.add_argument(
            "--breakpoint-cost-function",
            help="Distance used to compute the cost matrix (default: %(default)s). Options: "
            + ", ".join(d.BP_COST_FUNCTION),
            required=False,
            default=d.EUCLIDIAN,
            type=str,
        )
        optional_args_bp.add_argument(
            "--breakpoint-cost-function-threshold",
            help="Distance used to compute the cost matrix (default: %(default)s)",
            required=False,
            default=d.EUCLIDIAN_THRESHOLD,
            type=int,
        )

        optional_args_bp.add_argument(
            "--breakpoint-selection-distance",
            help="Distance for which two breakpoint-pairs are selected for matched candidates (default: %(default)s). Options: "
            + ", ".join(d.BP_CANDIDATES_SELECTION),
            required=False,
            default=d.MANHATTAN,
            type=str,
        )
        optional_args_bp.add_argument(
            "--breakpoint-selection-threshold",
            help="Threshold for which two breakpoint-pairs are considered matched candidates (default: %(default)s)",
            required=False,
            default=d.MANHATTAN_THRESHOLD,
            type=int,
        )
        # optional_args_bp.add_argument('--match-breakpoints-nonlinear',
        # 					help='Match breakpoints-pairs nonlinear, i.e. allow for high variability on one side if the other side matches good  (default: %(default)s). ' +
        # 					'Options: ' + ', '.join(d.BP_OPTIONS_NONLINEAR),
        # 					required=False, default=None, type=str)

        args = parser.parse_args()

        # 1. Get all configurations
        dict_configs = config(args)

        # 2. Create outdir if not exist
        if args.outdir is None:
            args.outdir = os.path.join(
                os.getcwd(), "output_" + str(datetime.datetime.now())
            )
            print(
                "Warning: outdir not specified. The files will be saved under ",
                args.outdir,
            )
        os.makedirs(args.outdir, exist_ok=True)

        # 3. Compare cycles
        if args.report and args.report == True:
            args.plot = True
        if args.plot is None:
            args.plot = False
        if args.report is None:
            args.report = False

        print("Structure A:", args.first_structure)
        print("Structure B:", args.second_structure)
        print("Outdir:", args.outdir)
        compare.compare_cycles(
            args.first_structure,
            args.second_structure,
            args.outdir,
            dict_configs,
            plot=args.plot,
            plot_report=args.report,
            min_cn=args.min_cn,
        )

    except Exception:
        traceback.print_exc()


if __name__ == "__main__":
    # pprint(Configs.DICT_DIST)
    main()
