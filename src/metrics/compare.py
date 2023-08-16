"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 12:05 PM 08/16/23

Compare cycle profiles.
"""

import os
import json

from metrics.features import *
from metrics.breakpoints import *
from metrics.distances import *
from utils import viz
from utils.utils import HEADER as ht
from utils.utils import DDT as ddt
from utils.utils import NpEncoder

import warnings
warnings.filterwarnings("ignore")

def remove_ref2methods(dict_configs):
	"""
	Remove from dictionary all values which are pointers to methods
	"""
	for key in dict_configs[ht.CONFIGS]:
		if key in [ht.BREAKPOINT_DISTANCE, ht.GENOMIC_FOOTPRINT]:
			for d in dict_configs[ht.CONFIGS][key]:
				del dict_configs[ht.CONFIGS][key][d][ht.DEFINITION]
	return dict_configs


def get_total_cost(dict_metrics):
	"""
	Get total cost
	"""
	total_cost = 0
	total_cost_description = {}
	for key in dict_metrics[ht.CONFIGS]:
		if key in [ht.BREAKPOINT_DISTANCE, ht.GENOMIC_FOOTPRINT]:
			for d in dict_metrics[ht.CONFIGS][key]:
				# distance is enabled
				if dict_metrics[ht.CONFIGS][key][d][ht.ENABLE]:
					weight = dict_metrics[ht.CONFIGS][key][d][ht.WEIGHT]
					val = dict_metrics[ht.DISTANCES][d]
					total_cost += weight * val
					total_cost_description[d] = weight

	return total_cost, total_cost_description


def compare_cycles(t_file: str, r_file: str, outdir: str, dict_configs: dict):
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
	dict_metrics = {}
	dict_metrics[ht.CONFIGS] = dict_configs
	dict_metrics[ht.DISTANCES] = {}

	# 1. Genome binning based on the breakpoints union

	# as data.frame objects
	df_t, df_r = read_input(t_file, r_file)
	bins = bin_genome(df_t, df_r, margin_size=0)

	# as pyranges objects
	df_t_pr, df_r_pr, df_bins_pr = read_pyranges(df_t, df_r, bins)

	# 2. Compute hamming distance and others
	h = get_hamming_score(df_t_pr, df_r_pr, df_bins_pr)
	h_norm = get_hamming_score_norm(df_t_pr, df_r_pr, df_bins_pr)

	dict_metrics[ht.DISTANCES][ddt.HAMMING] = round(h,2)
	dict_metrics[ht.DISTANCES][ddt.HAMMING_NORM] = round(h_norm,2)

	# 3. Compute copy-number similarity
	cv_profile_t = get_feature_cn(df_t, bins)
	cv_profile_r = get_feature_cn(df_r, bins)

	cv_distance = get_cosine_distance_cn(cv_profile_t, cv_profile_r)
	dict_metrics[ht.DISTANCES][ddt.COSINE_DISTANCE] = round(cv_distance,2)

	# 4. Breakpoint matching
	jc, matches, breakpoint_matches = compute_breakpoint_distance(df_t, df_r, distance=ddt.RELATIVE_METRIC)
	dict_metrics[ht.DISTANCES][ddt.JACCARD_DISTANCE] = jc

	# 5. Penalize for cycles and fragments multiplicity (if the tool decompose in one or more cycles)
	overlap_fragments_distance = get_overlap_fragments_weighted(df_t_pr, df_r_pr, df_bins_pr)
	overlap_cycles_distance = get_overlap_cycles_weighted(df_t_pr, df_r_pr, df_bins_pr)
	dict_metrics[ht.DISTANCES][ddt.FRAGMENTS_DISTANCE] = round(overlap_fragments_distance,2)
	dict_metrics[ht.DISTANCES][ddt.CYCLES_DISTANCE] = round(overlap_cycles_distance,2)

	# 6. Stoichiometry (compute distance of transforming one permutation in the other)

	# 7. Merge results
	total_cost, total_cost_description = get_total_cost(dict_metrics)
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST] = round(total_cost,2)
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST_DESCRIPTION] = total_cost_description

	# 8. Output
	with open(os.path.join(outdir, 'metrics.json'), 'w', encoding='utf-8') as f:
		final_dict = remove_ref2methods(dict_metrics)
		json.dump(final_dict, f, ensure_ascii=False, indent=4, cls=NpEncoder)

	# plot total cost
	viz.draw_total_cost(dict_metrics, os.path.join(outdir, 'total_cost.png'))

	# save coverage profile
	cv_profile_r.to_csv(os.path.join(outdir, 'e2_coverage_profile.txt'), header=True, index=False, sep="\t")
	cv_profile_t.to_csv(os.path.join(outdir, 'e1_coverage_profile.txt'), header=True, index=False, sep="\t")
