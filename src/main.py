"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 10:30 AM 07/28/23

AmpliconComparison entry point
"""

import argparse
import logging
import os
import traceback
import sys
from pprint import pprint

import metrics
from metrics import dict_distance_paired_breakpoints_thresholds as ddt
from configs import Configs
from utils import DDT as d
from utils import HEADER as h

def config(args):
	"""
	Configure comparison based on user spefication
	"""

	dict = Configs.DICT_DIST

	# genomic footprint distance settings
	dict[h.GENOMIC_FOOTPRINT][d.HAMMING_NORM][h.ENABLE] = args.hamming_distance_norm
	dict[h.GENOMIC_FOOTPRINT][d.COSINE_DISTANCE][h.ENABLE] = args.cosine_distance
	dict[h.GENOMIC_FOOTPRINT][d.FRAGMENTS_DISTANCE][h.ENABLE] = args.fragments_overlap_norm
	dict[h.GENOMIC_FOOTPRINT][d.CYCLES_DISTANCE][h.ENABLE] = args.cycles_overlap_norm

	# breakpoint matching distance settings
	dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = args.euclidian_distance
	dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.THRESHOLD] = args.euclidian_distance_threshold
	dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = args.relative_distance
	dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.THRESHOLD] = args.relative_distance_threshold

	return dict


def main():
	try:
		parser = argparse.ArgumentParser(prog="AmpliconComparison",
										 description='AmpliconComparison - compare cycle sets',
										 add_help=True)
		parser.add_argument('-a', '--first-structure', help='First structure (bed format)', required=True)
		parser.add_argument('-b', '--second-structure', help='Second structure (bed format)', required=True)
		parser.add_argument('-d', '--outdir', help='Output directory', required=True)
		parser.add_argument('--cos-similarity', help='Cosine similarity between coverage tracks (default: %(default)s)',
							required=False, default=True, type=bool)
		parser.add_argument('--hamming-distance-norm',
							help='Hamming distance between genomic footprint. Recommended when no copy-number information available (default: %(default)s)',
							required=False, default=True, type=bool)
		parser.add_argument('--cosine-distance',
						   help="Cosine distance between the coverage profile.",
							required=False, default=True, type=bool)
		parser.add_argument('--fragments-overlap-norm',
							help="Quantify the distance between fragments.",
							required=False, default=True, type=bool)
		parser.add_argument('--cycles-overlap-norm',
							help="Quantify the distance between cycles.",
							required=False, default=True, type=bool)
		parser.add_argument('--euclidian-distance',
							help='Use euclidian distance for breakpoint-pair matching (default: %(default)s)',
							required=False, default=False, type=bool)
		parser.add_argument('--euclidian-distance-threshold',
							help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
							required=False, default=3000, type=float)
		parser.add_argument('--relative-distance',
							help='Relative distance score for breakpoint matching (default: %(default)s)',
							required=False, default=False, type=bool)
		parser.add_argument('--relative-distance-threshold',
							help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
							required=False, default=0.3, type=float)
		# parser.add_argument('--unmatched-distance',
		# 					help='Max distance for which two breakpoint-pairs are considered unmatched (default: %(default)s)',
		# 					required=False, default=ddt["unmatched"], type=float)
		args = parser.parse_args()

		# 1. Get all configurations
		dict_configs = config(args)
		# 2. Create outdir if not exist
		os.makedirs(args.outdir, exist_ok=True)
		# 3. Compare cycles
		metrics.compare_cycles(args.first_structure,
							   args.second_structure,
							   args.outdir,
							   dict_configs)


	except Exception:
		traceback.print_exc()


if __name__ == '__main__':
	# pprint(Configs.DICT_DIST)
	main()
