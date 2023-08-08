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

import metrics
from metrics import dict_distance_paired_breakpoints_thresholds as ddt


def config(args):
	"""
	Configure comparison based on user spefication
	"""

	dict = {"cosine-similarity": args.cos_similarity,
			"hamming-distance": args.hamming_distance,
			"euclidian-distance": args.euclidian_distance,
			"euclidian-distance-threshold": args.euclidian_distance_threshold,
			"relative-distance": args.relative_distance,
			"relative-distance-threshold": args.relative_distance_threshold
			}

	return dict


def main():
	try:
		parser = argparse.ArgumentParser(prog="AmpliconComparison",
										 description='AmpliconComparison - compare cycle sets')
		parser.add_argument('-a', '--first-structure', help='First structure (bed format)', required=True)
		parser.add_argument('-b', '--second-structure', help='Second structure (bed format)', required=True)
		parser.add_argument('-d', '--outdir', help='Output directory', required=True)
		parser.add_argument('--cos-similarity', help='Cosine similarity between coverage tracks (default: %(default)s)',
							required=False, default=True, type=bool)
		parser.add_argument('--hamming-distance',
							help='Hamming distance between genomic footprint. Recommended when no copy-number information available (default: %(default)s)',
							required=False, default=True, type=bool)
		parser.add_argument('--euclidian-distance',
							help='Use euclidian distance for breakpoint-pair matching (default: %(default)s)',
							required=False,
							default=True, type=bool)
		parser.add_argument('--euclidian-distance-threshold',
							help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
							required=False,
							default=3000, type=float)
		parser.add_argument('--relative-distance',
							help='Relative distance score for breakpoint matching (default: %(default)s)',
							required=False,
							default=True, type=bool)
		parser.add_argument('--relative-distance-threshold',
							help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
							required=False,
							default=0.3, type=float)
		parser.add_argument('--unmatched-distance',
							help='Max distance for which two breakpoint-pairs are considered unmatched (default: %(default)s)',
							required=False,
							default=ddt["unmatched"], type=float)

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
	main()
