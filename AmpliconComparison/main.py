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

from AmpliconComparison.configs import Configs
from AmpliconComparison.utils.utils import DDT as d
from AmpliconComparison.utils.utils import HEADER as h
from AmpliconComparison.metrics import compare

# from AmpliconComparison.utils.utils import DDT as d
# from AmpliconComparison.utils.utils import HEADER as h

def config(args):
	"""
	Configure comparison based on user spefication
	"""

	dict = Configs.DICT_DIST

	# genomic footprint distance settings
	dict[h.GENOMIC_FOOTPRINT][d.HAMMING_NORM][h.ENABLE] = args.cn_hamming_dist
	dict[h.GENOMIC_FOOTPRINT][d.COSINE_DISTANCE][h.ENABLE] = args.cn_cosine_dist
	dict[h.GENOMIC_FOOTPRINT][d.COPYNUMBER_JC][h.ENABLE] = args.cn_jc_dist
	dict[h.GENOMIC_FOOTPRINT][d.FRAGMENTS_DISTANCE][h.ENABLE] = args.fragments_dist
	dict[h.GENOMIC_FOOTPRINT][d.CYCLES_DISTANCE][h.ENABLE] = args.cycles_dist

	# breakpoint matching distance settings
	# dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = args.euclidian_distance
	# dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.THRESHOLD] = args.euclidian_distance_threshold
	# dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = args.relative_distance
	# dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.THRESHOLD] = args.relative_distance_threshold

	dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE] = False
	dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.THRESHOLD] = 1000
	dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = False
	dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.THRESHOLD] = 0.3
	if args.breakpoint_dist:
		dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.EUCLIDIAN

		# if dict[h.BREAKPOINT_DISTANCE][d.EUCLIDIAN][h.ENABLE]:
		# 	dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.EUCLIDIAN
		# else:
		# 	dict[h.BREAKPOINT_DISTANCE][h.DEFAULT] = d.RELATIVE_METRIC
		# 	dict[h.BREAKPOINT_DISTANCE][d.RELATIVE_METRIC][h.ENABLE] = True

	return dict


def main():
	try:
		parser = argparse.ArgumentParser(prog="AmpliconComparison",
										 description='AmpliconComparison - compare cycle sets',
										 add_help=True)
		parser.add_argument('-a', '--first-structure', help='First structure (bed format)', required=True)
		parser.add_argument('-b', '--second-structure', help='Second structure (bed format)', required=True)
		parser.add_argument('-d', '--outdir', help='Output directory', required=True)
		parser.add_argument('-p', '--plot', help='Plot coverage profiles', default=True, type=bool)
		parser.add_argument('-r', '--report', help='Generate report (this will set \'plot\' also on True)', default=True, type=bool)
		parser.add_argument('--cn-hamming-dist',
							help='Hamming distance between genomic footprint. Recommended when no copy-number information available (default: %(default)s)',
							required=False, default=True, type=bool)
		parser.add_argument('--cn-cosine-dist',
						   help="Cosine distance between the coverage profile.",
							required=False, default=True, type=bool)
		parser.add_argument('--cn-jc-dist',
							help="Min-max distance / Jaccard distance between the coverage profile.",
							required=False, default=True, type=bool)
		parser.add_argument('--fragments-dist',
							help="Quantify the distance between fragments.",
							required=False, default=True, type=bool)
		parser.add_argument('--cycles-dist',
							help="Quantify the distance between cycles.",
							required=False, default=True, type=bool)
		parser.add_argument('--breakpoint-dist',
							help="Quantify the distance between cycles.",
							required=False, default=True, type=bool)
		# parser.add_argument('--euclidian-distance',
		# 					help='Use euclidian distance for breakpoint-pair matching (default: %(default)s)',
		# 					required=False, default=True, type=bool)
		# parser.add_argument('--euclidian-distance-threshold',
		# 					help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
		# 					required=False, default=1000, type=float)
		# parser.add_argument('--relative-distance',
		# 					help='Relative distance score for breakpoint matching (default: %(default)s)',
		# 					required=False, default=False, type=bool)
		# parser.add_argument('--relative-distance-threshold',
		# 					help='Distance threshold to accept two breakpoint-pairs as matched  (default: %(default)s)',
		# 					required=False, default=0.3, type=float)
		# parser.add_argument('--unmatched-distance',
		# 					help='Max distance for which two breakpoint-pairs are considered unmatched (default: %(default)s)',
		# 					required=False, default=ddt["unmatched"], type=float)
		args = parser.parse_args()

		# 1. Get all configurations
		dict_configs = config(args)

		# 2. Create outdir if not exist
		if args.outdir is None:
			args.outdir = os.path.join(os.getcwd(),"output_" + str(datetime.datetime.now()))
			print("Warning: outdir not specified. The files will be saved under ", args.outdir)
		os.makedirs(args.outdir, exist_ok=True)

		# 3. Compare cycles
		if args.report:
			args.plot = True
		compare.compare_cycles(args.first_structure,
							   args.second_structure,
							   args.outdir,
							   dict_configs,
							   plot=args.plot,
							   plot_report=args.report)


	except Exception:
		traceback.print_exc()


if __name__ == '__main__':
	# pprint(Configs.DICT_DIST)
	main()
