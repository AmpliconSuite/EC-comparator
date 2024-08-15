"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 12:05 PM 08/16/23

Compare cycle profiles.
"""

import os
import json
import pprint
import matplotlib.pyplot as plt

from src.metrics.features import *
from src.metrics.breakpoints import *
from src.metrics.distances import *
from src.utils import viz
from src.utils.utils import HEADER as ht
from src.utils.utils import DDT as ddt
from src.utils.utils import OUTFILES as o
from src.utils.utils import NpEncoder
# from src.utils import report


import warnings
warnings.filterwarnings("ignore")

def remove_ref2methods(dict_configs):
	"""
	Remove from dictionary all values which are pointers to methods
	"""
	for key in dict_configs[ht.CONFIGS]:
		if key in [ht.BREAKPOINT_DISTANCE, ht.GENOMIC_FOOTPRINT]:
			for d in dict_configs[ht.CONFIGS][key]:
				if ht.DEFINITION in dict_configs[ht.CONFIGS][key][d]:
					del dict_configs[ht.CONFIGS][key][d][ht.DEFINITION]
	return dict_configs


def get_total_cost(dict_metrics):
	"""
	Get total cost
	"""

	total_cost = 0
	total_cost_description = {}
	for key in dict_metrics[ht.CONFIGS]:
		if key == ht.GENOMIC_FOOTPRINT:
			for d in dict_metrics[ht.CONFIGS][key]:
				# distance is enabled
				if dict_metrics[ht.CONFIGS][key][d][ht.ENABLE]:
					weight = dict_metrics[ht.CONFIGS][key][d][ht.WEIGHT]
					val = dict_metrics[ht.DISTANCES][d]
					total_cost += weight * val
					total_cost_description[d] = weight
		elif key == ht.BREAKPOINT_DISTANCE:
			def_dist = dict_metrics[ht.CONFIGS][key][ht.DEFAULT]
			weight = dict_metrics[ht.CONFIGS][key][def_dist][ddt.JACCARD_DISTANCE][ht.WEIGHT]
			val = dict_metrics[ht.DISTANCES][ddt.JACCARD_DISTANCE]
			total_cost += weight * val
			total_cost_description[ddt.JACCARD_DISTANCE] = weight

	return total_cost, total_cost_description


def compare_cycles(t_file, r_file, outdir, dict_configs, plot=True):
	"""
	Entrypoint: compare the distance between two cycle sets.
	Args:
		t_file (str): First cycle set
		r_file (str): Second cycle set
		outdir (str): Output directory
		dict_configs (dict):
							hamming (bool): Include copy-number similarity for the total distance if set to True or Hamming distance is copy-number not available
							viz (bool): Enable visualization of the comparison

	Returns:
		Dictionary with different distances
	"""
	if outdir:
		os.makedirs(outdir, exist_ok=True)

	# pprint.pprint(dict_configs)

	dict_metrics = {}
	dict_metrics[ht.CONFIGS] = dict_configs
	default_breakpoint_distance = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ht.DEFAULT]
	default_breakpoint_distance_threshold = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][default_breakpoint_distance][ht.THRESHOLD]

	dict_metrics[ht.DISTANCES] = {}

	# 1. Genome binning based on the breakpoints union

	# as data.frame objects
	df_t, df_r = read_input(t_file, r_file)
	bins, chrlist = bin_genome(df_t, df_r, margin_size=0)
	chr_offsets = get_chromosome_offset(df_t, df_r)

	# as pyranges objects
	df_t_pr, df_r_pr, df_bins_pr = read_pyranges(df_t, df_r, bins)
	# as dataframe
	df_ = overlap_vector(df_t_pr, df_r_pr, df_bins_pr)

	# 2. Compute hamming distance and others
	df_, h = get_hamming_score(df_)
	h_norm = get_hamming_score_norm(df_)

	dict_metrics[ht.DISTANCES][ddt.HAMMING] = round(h,2)
	dict_metrics[ht.DISTANCES][ddt.HAMMING_NORM] = round(h_norm,2)

	# 3. Compute copy-number similarity
	cn_profile_t = get_feature_cn(df_t, df_)
	cn_profile_r = get_feature_cn(df_r, df_)

	cv_distance = get_cosine_distance_cn(cn_profile_t, cn_profile_r)
	dict_metrics[ht.DISTANCES][ddt.COSINE_DISTANCE] = round(cv_distance,2)
	cv_distance = get_jc_distance_cn(cn_profile_t, cn_profile_r)
	dict_metrics[ht.DISTANCES][ddt.COPYNUMBER_JC] = round(cv_distance, 2)

	# 4. Breakpoint matching
	br_t, br_r = get_breakpoints_pairs(df_t, df_r)
	jc, matches, breakpoint_matches = compute_breakpoint_distance(br_t, br_r,
																  distance=default_breakpoint_distance,
																  threshold=default_breakpoint_distance_threshold)
	dict_metrics[ht.DISTANCES][ddt.JACCARD_DISTANCE] = round(jc,2)
	# dict_metrics[ht.DISTANCES][ddt.JACCARD_DISTANCE] = 1

	# 5. Penalize for cycles and fragments multiplicity (if the tool decompose in one or more cycles)
	# overlap_fragments_distance = get_overlap_fragments_weighted(df_t_pr, df_r_pr, df_bins_pr)
	# overlap_cycles_distance = get_overlap_cycles_weighted(df_t_pr, df_r_pr, df_bins_pr)
	overlap_fragments_distance = get_overlap_fragments_weighted(df_)
	overlap_cycles_distance = get_overlap_cycles_weighted(df_t, df_r,  df_)
	dict_metrics[ht.DISTANCES][ddt.FRAGMENTS_DISTANCE] = round(overlap_fragments_distance,2)
	dict_metrics[ht.DISTANCES][ddt.CYCLES_DISTANCE] = round(overlap_cycles_distance,2)

	# 6. Stoichiometry (compute distance of transforming one permutation in the other)

	# 7. Merge results
	total_cost, total_cost_description = get_total_cost(dict_metrics)
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST] = round(total_cost,2)
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST_DESCRIPTION] = total_cost_description
	print("Total cost:",total_cost)

	# # 8. Output
	# if outdir:
	# 	with open(os.path.join(outdir, o.METRICS_JSON), 'w', encoding='utf-8') as f:
	# 		final_dict = remove_ref2methods(dict_metrics)
	# 		json.dump(final_dict, f, ensure_ascii=False, indent=4, cls=NpEncoder)
	#
	# # plot total cost
	# if plot:
	# 	outfile = os.path.join(outdir, o.TOTAL_COST_PNG) if outdir else None
	# 	viz.draw_total_cost(dict_metrics, outfile)
	# 	outfile = os.path.join(outdir, o.TOTAL_COST_TABLE) if outdir else None
	# 	viz.draw_total_cost_table(dict_metrics, outfile)
	#
	# # save and plot coverage profile
	# if outdir:
	# 	cn_profile_t.to_csv(os.path.join(outdir, o.COVERAGE_PROFILE_S1_TXT), header=True, index=False, sep="\t")
	# 	cn_profile_r.to_csv(os.path.join(outdir, o.COVERAGE_PROFILE_S2_TXT), header=True, index=False, sep="\t")
	# if plot:
	# 	outfile = os.path.join(outdir, o.COVERAGE_PROFILE_PNG) if outdir else None
	# 	viz.draw_cn(cn_profile_t, cn_profile_r, chrlist, outfile=outfile)
	#
	# # save breakpoints profile
	# if outdir:
	# 	br_t.to_csv(os.path.join(outdir, o.BREAKPOINTS_PROFILE_S1_TXT), header=True, index=False, sep="\t")
	# 	br_r.to_csv(os.path.join(outdir, o.BREAKPOINTS_PROFILE_S2_TXT), header=True, index=False, sep="\t")
	#
	# # plot merged coverage and breakpoint profile
	# max_coverage = max(np.max(cn_profile_t[ht.CN].tolist()), np.max(cn_profile_r[ht.CN].tolist()))
	# max_coverage = max_coverage + 0.5 * max_coverage
	# if plot:
	# 	outfile = os.path.join(outdir, o.COVERAGE_BREAKPOINTS_PROFILE) if outdir else None
	# 	viz.plot_combined(br_t, br_r, cn_profile_t, cn_profile_r, breakpoint_matches, chrlist, max_coverage, outfile=outfile)
	#
	# # create report
	# if outdir:
	# 	report.generate_report(
	# 		t_file,
	# 		r_file,
	# 		outdir=outdir,
	# 		total_cost=total_cost)