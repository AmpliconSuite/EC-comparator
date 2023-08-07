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


def main():
	try:
		parser = argparse.ArgumentParser(prog="AmpliconComparison",
										 description='AmpliconComparison - compare cycle sets')
		parser.add_argument('-a', '--first-structure', help='First structure (bed format)', required=True)
		parser.add_argument('-b', '--second-structure', help='Second structure (bed format)', required=True)
		parser.add_argument('-d', '--outdir', help='Output directory', required=True)
		parser.add_argument('--cos-similarity', help='Cosine similarity between coverage tracks (default: %(default)s',
							required=False, default=True, type=bool)
		parser.add_argument('--hamming-distance',
							help='Hamming distance between genomic footprint. Recommended when no copy-number information available (default: %(default)s',
							required=False, default=True, type=bool)
		parser.add_argument('--euclidian-distance',
							help='Use euclidian distance for breakpoint matching (default: %(default)s', required=False,
							default=True, type=bool)
		parser.add_argument('--relative-distance-score',
							help='Relative distance score for breakpoint matching (default: %(default)s',
							required=False,
							default=True, type=bool)
		args = parser.parse_args()
		print(args)
	except Exception:
		traceback.print_exc()


if __name__ == '__main__':
	main()
