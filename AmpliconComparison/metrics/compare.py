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

from AmpliconComparison.metrics.features import *
from AmpliconComparison.metrics.breakpoints import *
from AmpliconComparison.metrics.distances import *
from AmpliconComparison.utils import viz
from AmpliconComparison.utils.utils import HEADER as ht
from AmpliconComparison.utils.utils import DDT as ddt
from AmpliconComparison.utils.utils import OUTFILES as o
from AmpliconComparison.utils.utils import NpEncoder
from AmpliconComparison.utils import report
from AmpliconComparison.utils.utils import get_weight_distance, get_value_distance


import warnings
warnings.filterwarnings("ignore")

def remove_key(d, key_to_remove='definition'):
	"""
	Remove from dictionary all values which are pointers to methods
	"""
	if isinstance(d, dict):
		return {k: remove_key(v, key_to_remove) for k, v in d.items() if k != key_to_remove}
	elif isinstance(d, list):
		return [remove_key(v, key_to_remove) for v in d]
	else:
		return d



def get_total_cost(dict_metrics):
	"""
	Get total cost
	"""

	total_cost = 0
	total_cost_description = {}

	for key in dict_metrics[ht.CONFIGS]:

		if key == ht.GENOMIC_FOOTPRINT:
			for d in dict_metrics[ht.CONFIGS][key]:
				weight, _ = get_weight_distance(key,d,dict_metrics)
				val, _ = get_value_distance(key,d,dict_metrics)
				total_cost += weight * val
				total_cost_description[d] = weight

		elif key == ht.BREAKPOINT_DISTANCE:
			# take the weight and value of the default distance for breakpoint-matching
			weight, d = get_weight_distance(key, None, dict_metrics)
			val, d = get_value_distance(key, None, dict_metrics)
			total_cost += weight * val
			total_cost_description[d] = weight

	return total_cost, total_cost_description


def compare_cycles(t_file, r_file, outdir, dict_configs, plot=True, plot_report=True, min_cn=0, s1=None,s2=None):
	"""
	Entrypoint: compare the distance between two structures sets.
	Args:
		t_file (str): First structures set
		r_file (str): Second structures set
		outdir (str): Output directory
		dict_configs (dict):
							hamming (bool): Include copy-number similarity for the total distance if set to True or Hamming distance is copy-number not available
							viz (bool): Enable visualization of the comparison
       plot (bool): Generate coverage and breakpoints plot comparatively for first and second structure
       plot_report (bool): Generate a report including a summary plot of the distances and coverage and breakpoints plots
       min_cn (float): Set a minimal copy-number or coverage for the structures to be considered for comparison (default 0)
       s1 (str): Name of structure 1
       s2 (str): Name of structure 2

	Returns:
		Dictionary with different distances
	"""
	if outdir:
		os.makedirs(outdir, exist_ok=True)

	# pprint.pprint(dict_configs)
	# set names
	s1 = os.path.basename(t_file)
	s2 = os.path.basename(r_file)

	dict_metrics = {}
	dict_metrics[ht.CONFIGS] = dict_configs
	
	
	# umatching breakpoints 
	default_unmatching_distance = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ddt.UNMATCHED_DISTANCE]
	default_unmatching_threshold = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ddt.UNMATCHED_THRESHOLD]
	
	# compute cost matrix using distance d
	default_breakpoint_distance = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ht.DEFAULT]
	default_breakpoint_distance_threshold = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][default_breakpoint_distance][ht.THRESHOLD]
 
	# allow nonlinear matching of breakpoints
	# default_match_nonlinear = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ddt.MATCH_NONLINEAR]
 
	# how to compute the breakpoint matching distance
	default_breakpoint_distance_calculation = dict_metrics[ht.CONFIGS][ht.BREAKPOINT_DISTANCE][ht.BREAKPOINT_DISTANCE_CALCULATION]
	dict_metrics[ht.DISTANCES] = {}

	# 1. Genome binning based on the breakpoints union

	# as data.frame objects
	df_t, df_r = read_input(t_file, r_file, outdir, min_cn=min_cn)
	
	if df_t.shape[0] == 0:
		raise Exception("""EC-comparator stoppe because of empty file: {} is empty""".format(t_file))

	if df_r.shape[0] == 0:
		raise Exception("""EC-comparator stoppe because of empty file: {} is empty""".format(r_file))
	
	# print(df_t)
	# print(df_r)
 
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
																  distance_threshold=default_breakpoint_distance_threshold,
																  unmatched_dist=default_unmatching_distance,
																  unmatched_threshold=default_unmatching_threshold,
																  how=default_breakpoint_distance_calculation,
                  												#   match_nonlinear=default_match_nonlinear
                                )
	print("Breakpoint maching: JD:", jc)
	dict_metrics[ht.DISTANCES][ddt.JACCARD_DISTANCE] = round(jc,2)
 
	# 5. Penalize for cycles and fragments multiplicity (if the tool decompose in one or more cycles)
	# overlap_fragments_distance = get_overlap_fragments_weighted(df_t_pr, df_r_pr, df_bins_pr)
	# overlap_cycles_distance = get_overlap_cycles_weighted(df_t_pr, df_r_pr, df_bins_pr)
	overlap_fragments_distance = get_overlap_fragments_weighted(df_)
	overlap_cycles_distance = get_overlap_cycles_weighted(df_t, df_r,  df_)
	dict_metrics[ht.DISTANCES][ddt.FRAGMENTS_DISTANCE] = round(overlap_fragments_distance,2)
	dict_metrics[ht.DISTANCES][ddt.CYCLES_DISTANCE] = round(overlap_cycles_distance,2)

	# 6. Stoichiometry (compute distance of transforming one permutation in the other)
	# missing

	# 7. Merge results
	total_cost, total_cost_description = get_total_cost(dict_metrics)
	total_cost = round(total_cost,2)
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST] = total_cost
	dict_metrics[ht.DISTANCES][ddt.TOTAL_COST_DESCRIPTION] = total_cost_description
	print("Total cost:",total_cost)
	print("Plot:",plot)
	print("Report:",plot_report)

	# 8. Output
	if outdir:
		with open(os.path.join(outdir, o.METRICS_JSON), 'w', encoding='utf-8') as f:
			final_dict = remove_key(dict_metrics,key_to_remove='definition')
			json.dump(final_dict, f, ensure_ascii=False, indent=4, cls=NpEncoder)

		# save coverage profile
		cn_profile_t.to_csv(os.path.join(outdir, o.COVERAGE_PROFILE_S1_TXT), header=True, index=False, sep="\t")
		cn_profile_r.to_csv(os.path.join(outdir, o.COVERAGE_PROFILE_S2_TXT), header=True, index=False, sep="\t")

		# save breakpoint profile
		br_t.to_csv(os.path.join(outdir, o.BREAKPOINTS_PROFILE_S1_TXT), header=True, index=False, sep="\t")
		br_r.to_csv(os.path.join(outdir, o.BREAKPOINTS_PROFILE_S2_TXT), header=True, index=False, sep="\t")

		if plot:

			# plot coverage profile
			outfile = os.path.join(outdir, o.COVERAGE_PROFILE_PNG)
			viz.draw_cn(cn_profile_t, cn_profile_r, chrlist,
						outfile=outfile,s1=s1,s2=s2)


			# plot merged coverage and breakpoint profile
			max_coverage = max(np.max(cn_profile_t[ht.CN].tolist()), np.max(cn_profile_r[ht.CN].tolist()))
			max_coverage = max_coverage + 0.5 * max_coverage

			outfile = os.path.join(outdir, o.COVERAGE_BREAKPOINTS_PROFILE) if outdir else None
			viz.plot_combined(br_t, br_r, cn_profile_t, cn_profile_r, breakpoint_matches,
							  chrlist, max_coverage, outfile=outfile,s1=s1,s2=s2)

			# plot total cost
			outfile = os.path.join(outdir, o.TOTAL_COST_PNG)
			viz.draw_total_cost(dict_metrics, outfile)
			outfile = os.path.join(outdir, o.TOTAL_COST_TABLE)
			viz.draw_total_cost_table(dict_metrics, outfile)

			# create report
			if plot_report:
				report.generate_report(
					t_file,
					r_file,
					outdir=outdir,
					total_cost=total_cost)