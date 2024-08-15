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

from src.utils.utils import HEADER as ht
from src.utils.utils import DDT as ddt
from src.utils.utils import PROPS as p
from src.metrics.features import rename_columns

def get_hamming_score(df_):
	"""
	Compute hamming distance between the two reconstructions

	Args:
		df_(pd.DataFrame):

	"""
	df_[ddt.HAMMING] = df_.apply(lambda x: np.logical_xor(x.h1, x.h2), axis=1)
	df_[ddt.HAMMING] = df_[ddt.HAMMING].astype(int)
	df_[ddt.PROD] = df_[ddt.HAMMING] * (df_[ht.CEND] - df_[ht.CSTART])
	df_[ddt.LEN] = abs(df_[ht.CEND] - df_[ht.CSTART]) * df_[ht.BIN_ENABLED]

	return df_, df_[ddt.PROD].sum()


def get_hamming_score_norm(df_):
	"""
	Compute hamming distance between the two reconstructions

	Args:
		df_(pd.DataFrame)
	"""

	# df_['hamming'] = df_.apply(lambda x: np.logical_xor(x.h1, x.h2), axis=1)
	# df_['hamming'] = df_['hamming'].astype(int)
	#
	# # consider only bins used in the reconstruction of s1 or s2
	# df_["len"] = abs(df_["End"] - df_["Start"]) * df_["bins"]
	# df_["prod"] = df_["hamming"] * df_["len"]

	return (df_[ddt.PROD].sum()) / (df_[ddt.LEN].sum())


def get_overlap_fragments_weighted(df_):
	"""
	Compute the overlap distance between the fragment counts of two reconstructions which overlap every genomic bin.

	Args:
		df_ (pd.DataFrame): Binned summary reconstructions
	"""
	df_[ddt.OVERLAP_SCORE] = df_.apply(lambda x: abs(x.e1 - x.e2), axis=1)
	# important to disable unused genomic bins
	df_[ddt.LEN] = abs(df_[ht.CEND] - df_[ht.CSTART]) * df_[ht.BIN_ENABLED]

	# account for duplicated fragments
	df_[ddt.TOTALLEN] = df_.apply(lambda x: max(x.e1, x.e2) * x.len, axis=1)
	df_[ddt.PROD] = df_[ddt.OVERLAP_SCORE] * df_[ddt.LEN]

	return (df_[ddt.PROD].sum()) / (df_[ddt.TOTALLEN].sum())

def count_cycles(results, e1,bins,c, col):
	"""
	Args:
		results (pd.DataFrame):
		e1 (pd.DataFrame): Amplicon reconstructions
		bins (pd.DataFrame): Reference genomic bins
		c (str): Cycle_id
		col (str): Has value 's1' or 's2'
	"""

	# filter pyrange by circle name
	df_tmp = e1[e1[ht.CIRC_ID] == c]
	df_tmp_pr = pr.PyRanges(
		df_tmp.rename(columns=rename_columns(df_tmp.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))

	# compute overlap
	grs = {n: s for n, s in zip([ht.BINS, col], [bins, df_tmp_pr])}
	overlap = pr.count_overlaps(grs).as_df()
	# mark bin as covered by cycle c
	overlap[col] = overlap[col].apply(lambda x: 1 if x >= 1 else 0)

	# merge results
	results = pd.concat([results, overlap]).groupby(
		[ht.CCHR, ht.CSTART, ht.CEND], observed=True, as_index=False).sum()

	return results

def get_overlap_cycles_weighted(e1, e2, bins):
	"""
	Compute the overlap distance between the cycles counts of two reconstructions which overlap every genomic bin.

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DataFrame): Bins of the intervals union e1 and e2
	"""
	pr_bins = pr.PyRanges(bins[[ht.CCHR,ht.CSTART,ht.CEND]])
	c1 = e1[ht.CIRC_ID].drop_duplicates().tolist()
	c2 = e2[ht.CIRC_ID].drop_duplicates().tolist()

	overlap_parent1 = deepcopy(bins)
	overlap_parent1[ht.S1] = 0

	# overlap_parent1[ddt.LEN] = abs(overlap_parent1[ht.CEND] - overlap_parent1[ht.CSTART])
	for c in c1:
		overlap_parent1 = count_cycles(overlap_parent1, e1, pr_bins, c, ht.S1)

	overlap_parent2 = deepcopy(bins)
	overlap_parent2[ht.S2] = 0

	for c in c2:
		overlap_parent2 = count_cycles(overlap_parent2, e2, pr_bins, c, ht.S2)

	overlaps = pd.merge(overlap_parent1[[ht.CCHR, ht.CSTART, ht.CEND, ht.S1, ht.BIN_ENABLED]],
						overlap_parent2[[ht.CCHR, ht.CSTART, ht.CEND, ht.S2, ht.BIN_ENABLED]], how='inner')

	overlaps[ddt.OVERLAP_SCORE] = overlaps.apply(lambda x: abs(x[ht.S1] - x[ht.S2]), axis=1)
	overlaps[ddt.LEN] = abs(overlaps[ht.CEND] - overlaps[ht.CSTART]) * overlaps[ht.BIN_ENABLED]
	overlaps[ddt.TOTALLEN] = overlaps.apply(lambda x: max(x[ht.S1], x[ht.S2]) * x.len, axis=1)
	overlaps[ddt.PROD] = overlaps[ddt.OVERLAP_SCORE] * overlaps[ddt.LEN]

	return (overlaps[ddt.PROD].sum()) / (overlaps[ddt.TOTALLEN].sum())


def get_cosine_distance_cn(cn_profile1, cn_profile2):
	"""
	Get cosine distance between the two copy number profiles
	"""
	cn_profile1["estimated_cn_normalized"] = cn_profile1.apply(lambda x: x[ht.CN] * abs(x[ht.CEND] - x[ht.CSTART]),
															   axis=1)
	cn_profile2["estimated_cn_normalized"] = cn_profile2.apply(lambda x: x[ht.CN] * abs(x[ht.CEND] - x[ht.CSTART]),
															   axis=1)

	# filter out bins which are disabled
	cn_profile1_new = cn_profile1[cn_profile1[ht.BIN_ENABLED] == 1]
	cn_profile2_new = cn_profile2[cn_profile2[ht.BIN_ENABLED] == 1]

	a = np.array([cn_profile1_new["estimated_cn_normalized"].tolist()])
	b = np.array([cn_profile2_new["estimated_cn_normalized"].tolist()])

	return 1 - abs(cosine_similarity(a,b))[0][0]


def get_jc_distance_cn(cn_profile1, cn_profile2):
	"""
	Get jaccard index distance between the two copy number profiles
	"""
	cn_merge = pd.merge(cn_profile1, cn_profile2, on=[ht.CCHR, ht.CSTART, ht.CEND, ht.BIN_ENABLED])
	cn_merge["max"] = cn_merge.apply(lambda x: max(x["estimated_cn_normalized_x"], x["estimated_cn_normalized_y"]),
									 axis=1)
	# overlap
	cn_merge["min"] = cn_merge.apply(lambda x: min(x["estimated_cn_normalized_x"], x["estimated_cn_normalized_y"]),
									 axis=1)

	# filter out unused bins
	cn_merge = cn_merge[cn_merge[ht.BIN_ENABLED]==1]

	# JC = matches / all
	# 1 - JC means -> 0 highly similar, 1 - not similar
	return 1 - cn_merge["min"].sum() / cn_merge["max"].sum()


def euclidian_distance(cha, a, chb, b, chx, x, chy, y):
	if cha == chx and chb == chy:
		return math.sqrt((a - x) ** 2 + (b - y) ** 2)
	elif cha == chy and chb == chx:
		return math.sqrt((a - y) ** 2 + (b - x) ** 2)
	else:
		return np.nan


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
						  ddt.OVERLAP: get_overlap_fragments_weighted,
						  ddt.COSINE_DISTANCE: get_cosine_distance_cn,
						  ddt.CYCLES_DISTANCE: get_overlap_cycles_weighted,
						  ddt.FRAGMENTS_DISTANCE: get_overlap_fragments_weighted,
						  ddt.COPYNUMBER_JC: get_jc_distance_cn}

# distance between paired breapoints
dict_distance_paired_breakpoints = {
	ddt.EUCLIDIAN: euclidian_distance,
	# ddt.EUCLIDIAN_NORM_L1: euclidian_distance_norm_l1,
	# ddt.EUCLIDIAN_NORM_L2: euclidian_distance_norm_l2,
	# "auc_triangle": auc_triangle,
	# "auc_trapeze": auc_trapeze,
	# "match_angle": match_angle,
	# ddt.RELATIVE_METRIC: match_score
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

