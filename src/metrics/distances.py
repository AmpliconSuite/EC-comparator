"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 11:45 AM 08/16/23

Compute different distances for a set of features.
"""
import math
import pandas as pd
import numpy as np
from copy import deepcopy

from sklearn.metrics.pairwise import cosine_similarity
import sklearn.metrics.pairwise as skl
from sklearn import preprocessing

import pyranges as pr

from utils.utils import HEADER as ht
from utils.utils import DDT as ddt


def get_cos_similarity_cn(e1, e2):
	"""
	Args:
		e1 (pd.DataFrame): First reconstruction binned by union of all bins
		e2 (pd.DataFrame): Second reconstruction binned by union of all bins
	"""
	return cosine_similarity(e1["cn"].tolist(), e2["cn"].tolist())


def get_hamming_score(e1, e2, bins):
	"""
	Compute hamming distance between the two reconstructions

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DataFrame): Bins of the intervals union e1 and e2
	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()

	# everything > 1 set to 1
	overlaps['e1'] = overlaps['e1'].apply(lambda x: 1 if x >= 1 else 0)
	overlaps['e2'] = overlaps['e2'].apply(lambda x: 1 if x >= 1 else 0)

	overlaps['hamming'] = overlaps.apply(lambda x: np.logical_xor(x.e1, x.e2), axis=1)
	overlaps['hamming'] = overlaps['hamming'].astype(int)

	overlaps["prod"] = overlaps["hamming"] * (overlaps["End"] - overlaps["Start"])
	return overlaps["prod"].sum()


def get_hamming_score_norm(e1, e2, bins):
	"""
	Compute hamming distance between the two reconstructions

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DatamFrame): Bins of the intervals union e1 and e2
	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()

	# everything > 1 set to 1
	overlaps['e1'] = overlaps['e1'].apply(lambda x: 1 if x >= 1 else 0)
	overlaps['e2'] = overlaps['e2'].apply(lambda x: 1 if x >= 1 else 0)

	overlaps['hamming'] = overlaps.apply(lambda x: np.logical_xor(x.e1, x.e2), axis=1)
	overlaps['hamming'] = overlaps['hamming'].astype(int)

	overlaps["len"] = abs(overlaps["End"] - overlaps["Start"])
	overlaps["prod"] = overlaps["hamming"] * overlaps["len"]

	return (overlaps["prod"].sum()) / (overlaps["len"].sum())


def get_overlap_fragments_weighted(e1, e2, bins):
	"""
	Compute the overlap distance between the fragment counts of two reconstructions which overlap every genomic bin.

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DatamFrame): Bins of the intervals union e1 and e2
	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()
	overlaps["overlapping_score"] = overlaps.apply(lambda x: abs(x.e1 - x.e2), axis=1)
	overlaps["len"] = abs(overlaps["End"] - overlaps["Start"])
	overlaps["prod"] = overlaps["overlapping_score"] * overlaps["len"]

	return (overlaps["prod"].sum()) / (overlaps["len"].sum())


def get_overlap_cycles_weighted(e1, e2, bins):
	"""
	Compute the overlap distance between the cycles counts of two reconstructions which overlap every genomic bin.

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DatamFrame): Bins of the intervals union e1 and e2
	"""
	c1 = e1.as_df()[ht.CIRC_ID].drop_duplicates().tolist()
	c2 = e2.as_df()[ht.CIRC_ID].drop_duplicates().tolist()

	overlap_parent1 = deepcopy(bins.as_df())
	overlap_parent1["e1"] = 0
	for c in c1:
		# filter pyrange by circle name
		df_tmp = e1.as_df()
		df_tmp = df_tmp[df_tmp[ht.CIRC_ID] == c]
		df_tmp_pr = pr.PyRanges(df_tmp)

		# compute overlap
		grs = {n: s for n, s in zip(["bins", "e1"], [bins, df_tmp_pr])}
		overlap = pr.count_overlaps(grs).as_df()
		overlap['e1'] = overlap['e1'].apply(lambda x: 1 if x >= 1 else 0)

		# merge results
		overlap_parent1 = pd.concat([overlap, overlap_parent1]).groupby(
			["Chromosome", "Start", "End", "bins"]).sum().reset_index()

	overlap_parent2 = deepcopy(bins.as_df())
	overlap_parent2["e2"] = 0
	for c in c2:
		# filter pyrange by circle name
		df_tmp = e2.as_df()
		df_tmp = df_tmp[df_tmp[ht.CIRC_ID] == c]
		df_tmp_pr = pr.PyRanges(df_tmp)

		# compute overlap
		grs = {n: s for n, s in zip(["bins", "e2"], [bins, df_tmp_pr])}
		overlap = pr.count_overlaps(grs).as_df()
		overlap['e2'] = overlap['e2'].apply(lambda x: 1 if x >= 1 else 0)

		# merge results
		overlap_parent2 = pd.concat([overlap, overlap_parent2]).groupby(
			["Chromosome", "Start", "End", "bins"]).sum().reset_index()

	overlaps = pd.merge(overlap_parent1[["Chromosome", "Start", "End", "e1"]],
						overlap_parent2[["Chromosome", "Start", "End", "e2"]], how='inner')
	overlaps["overlapping_score"] = overlaps.apply(lambda x: abs(x.e1 - x.e2), axis=1)
	overlaps["len"] = abs(overlaps["End"] - overlaps["Start"])
	overlaps["prod"] = overlaps["overlapping_score"] * overlaps["len"]

	return (overlaps["prod"].sum()) / (overlaps["len"].sum())


def get_cosine_distance_cn(cn_profile1, cn_profile2):
	"""
	Get cosine distance between the two copy number profiles
	"""
	return 1 - abs(cosine_similarity(np.array([cn_profile1[ht.CN].tolist()]),
									 np.array([cn_profile2[ht.CN].tolist()]))[0][0])


def euclidian_distance(a, b, x, y):
	return math.sqrt((a - x) ** 2 + (b - y) ** 2)


def euclidian_distance_norm_l2(a, b, x, y):
	v1 = np.array([[a, b]])
	v2 = np.array([[x, y]])
	v1_norm = preprocessing.normalize(v1, norm='l2')
	v2_norm = preprocessing.normalize(v2, norm='l2')

	return skl.euclidean_distances(v1_norm, v2_norm)


def euclidian_distance_norm_l1(a, b, x, y):
	v1 = np.array([[a, b]])
	v2 = np.array([[x, y]])
	v1_norm = preprocessing.normalize(v1, norm='l1')
	v2_norm = preprocessing.normalize(v2, norm='l1')
	return skl.euclidean_distances(v1_norm, v2_norm)


def auc_triangle(a, b, x, y):
	return abs(a - x) * abs(b - y) / 2


def auc_trapeze(a, b, x, y):
	return (b + y) * abs(a - x) / 2


def match_angle(a, b, x, y):
	alpha = math.atan(abs(x - y) / abs(a - b))
	alpha_45 = math.pi / 4
	# return 1 - alpha / alpha_45
	return alpha / alpha_45

def match_score(a, b, x, y):
	"""
	Compute the similarity between relative distances of the breakpoints.

	Arguments:
		a (float):
		b (float):
		x (float):
		y (float:
	Returns:
		score
	"""

	v1 = np.array([[a - x,
					a - b,
					a - y]])
	v2 = np.array([[x - a,
					x - b,
					x - y]])
	v3 = np.array([[b - a,
					b - x,
					b - y]])
	v4 = np.array([[y - a,
					y - x,
					y - b]])

	cos1 = cosine_similarity(v1, v2)
	cos2 = cosine_similarity(v3, v4)

	# return (cos1 + cos2) / 2
	return 1 - abs((cos1 + cos2) / 2)


# distance for the copy-number profile
dict_distance_function = {ddt.HAMMING: get_hamming_score,
						  ddt.HAMMING_NORM: get_hamming_score_norm,
						  ddt.OVERLAP: get_overlap_fragments_weighted}

# distance between paired breapoints
dict_distance_paired_breakpoints = {
	ddt.EUCLIDIAN: euclidian_distance,
	# ddt.EUCLIDIAN_NORM_L1: euclidian_distance_norm_l1,
	ddt.EUCLIDIAN_NORM_L2: euclidian_distance_norm_l2,
	# "auc_triangle": auc_triangle,
	# "auc_trapeze": auc_trapeze,
	# "match_angle": match_angle,
	ddt.RELATIVE_METRIC: match_score
}

dict_distance_paired_breakpoints_thresholds = {
	ddt.EUCLIDIAN: ddt.EUCLIDIAN_THRESHOLD,
	ddt.EUCLIDIAN_NORM_L1: ddt.EUCLIDIAN_THRESHOLD_L1,
	ddt.EUCLIDIAN_NORM_L2: ddt.EUCLIDIAN_THRESHOLD_L2,
	ddt.AUC_TRIANGLE: ddt.AUC_TRIANGLE_THRESHOLD,
	ddt.AUC_TRAPEZE: ddt.AUC_TRAPEZE_THRESHOLD,
	ddt.RELATIVE_METRIC: ddt.RELATIVE_METRIC_THRESHOLD,
	ddt.MATCH_ANGLE: ddt.MATCH_ANGLE_THRESHOLD,
	ddt.UNMATCHED: ddt.UNMATCHED_THRESHOLD
}
