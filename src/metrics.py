import os
import collections
import math
import json
import pandas as pd
import numpy as np
from copy import deepcopy
from pprint import pprint

from sklearn.metrics.pairwise import cosine_similarity
import sklearn.metrics.pairwise as skl
from sklearn import preprocessing

from pybedtools import BedTool
import pyranges as pr
import networkx as nx

from utils import HEADER as ht
from utils import PROPS as p
from utils import DDT as ddt
from utils import NpEncoder


def rename_columns(df_cols, dict_mapping_cols):
	"""
	Try to rename columns so that it matched the naming
	"""
	dict_newcols = {}
	for c in df_cols:
		if c in dict_mapping_cols:
			dict_newcols[c] = dict_mapping_cols[c]
		else:
			dict_newcols[c] = c
	return dict_newcols


def read_input(t_file, r_file):
	"""
	Read input true and reconstruct file
	"""
	dtype = {'#chr': 'str',
			 'chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'cycle_id': 'str',
			 'weight': 'double',
			 'score': 'double',
			 'estimated_cn': 'double',
			 'strand': 'str',
			 'orientation': str}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	# rename columns and check if all columns available
	t_collection.rename(columns=rename_columns(t_collection.columns.tolist(), ht.DICT_HEADER), errors="raise",
						inplace=True)
	r_collection.rename(columns=rename_columns(r_collection.columns.tolist(), ht.DICT_HEADER), errors="raise",
						inplace=True)

	# reorder columns
	lorder = ht.HEADER_SORTED
	for h in t_collection.columns.tolist():
		if h not in lorder:
			lorder.append(h)

	return t_collection[lorder], r_collection[lorder]


def bin_genome(t_collection, r_collection, margin_size=10000):
	"""
	Bin the intervals by the breakpoints union.
	Warning: Setting a margin size can effect the output of different distances
	"""
	df_bins = pd.DataFrame(np.concatenate((t_collection[[ht.CHR, ht.START]].values,
										   t_collection[[ht.CHR, ht.END]].values,
										   r_collection[[ht.CHR, ht.START]].values,
										   r_collection[[ht.CHR, ht.END]].values),
										  axis=0))

	df_bins.columns = [ht.CHR, ht.START]
	# return df_bins

	# get max and min for each chromosome
	df_bins_gr = df_bins.groupby([ht.CHR]).agg({ht.START: [np.min, np.max]}).reset_index()
	df_bins_gr.columns = [ht.CHR, "start_amin", "start_amax"]
	df_bins_gr["start_amin"] = df_bins_gr.apply(lambda x: max(0, x["start_amin"] - margin_size), axis=1)
	df_bins_gr["start_amax"] = df_bins_gr.apply(lambda x: x["start_amax"] + margin_size, axis=1)

	df_bins = pd.concat([df_bins, df_bins_gr[[ht.CHR, "start_amin"]].rename(columns={"start_amin": ht.START})],
						ignore_index=True)
	df_bins = pd.concat([df_bins, df_bins_gr[[ht.CHR, "start_amax"]].rename(columns={"start_amax": ht.START})],
						ignore_index=True)
	df_bins = df_bins.sort_values(by=[ht.CHR, ht.START]).drop_duplicates()

	# rotate with 1 up the start column
	df_bins_suffix = df_bins.tail(-1)
	df_bins_suffix = df_bins_suffix.append((df_bins.head(1)), ignore_index=True)

	df_bins.reset_index(drop=True, inplace=True)
	df_bins_suffix.reset_index(drop=True, inplace=True)
	df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
	df_bins.columns = [ht.CHR, ht.START, "#chr2", ht.END]
	df_bins[ht.LEN] = df_bins[ht.END].astype(int) - df_bins[ht.START].astype(int)

	# keep only rows with same chr and positive distance
	df_bins = df_bins[(df_bins[ht.CHR] == df_bins["#chr2"]) & (df_bins[ht.LEN] > 0)]

	return df_bins[[ht.CHR, ht.START, ht.END, ht.LEN]]


# def bin_genome(t_collection, r_collection, margin_size=10000):
# 	"""
# 	Bin the intervals by the breakpoints union.
# 	"""
# 	df_bins = pd.DataFrame(np.concatenate((t_collection[[ht.CHR, ht.START]].values,
# 										   t_collection[[ht.CHR, ht.END]].values,
# 										   r_collection[[ht.CHR, ht.START]].values,
# 										   r_collection[[ht.CHR, ht.END]].values),
# 										  axis=0))
# 	df_bins.columns = [ht.CHR, ht.START]
# 	df_bins = df_bins.sort_values(by=[ht.CHR, ht.START]).drop_duplicates()
#
# 	# rotate with 1 up the start column
# 	df_bins_suffix = df_bins.tail(-1)
# 	df_bins_suffix = df_bins_suffix.append(df_bins.head(1), ignore_index=True)
#
# 	df_bins.reset_index(drop=True, inplace=True)
# 	df_bins_suffix.reset_index(drop=True, inplace=True)
# 	df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
# 	df_bins.columns = [ht.CHR, ht.START, "#chr2", ht.END]
# 	df_bins[ht.LEN] = df_bins[ht.END].astype(int) - df_bins[ht.START].astype(int)
#
# 	# keep only rows with same chr and positive distance
# 	df_bins = df_bins[(df_bins[ht.CHR] == df_bins["#chr2"]) & (df_bins[ht.LEN] > 0)]
#
# 	return df_bins[[ht.CHR, ht.START, ht.END, ht.LEN]]

def read_pyranges(t_file, r_file, bins):
	"""
	Load data as pyranges
	"""

	df_t_pr = pr.read_bed(t_file, as_df=False, nrows=None)
	df_r_pr = pr.read_bed(r_file, as_df=False, nrows=None)

	# df_t_pr = deepcopy(df_t)
	# df_r_pr = deepcopy(df_r)
	#
	# # rename columns
	# df_t_pr.rename(columns=rename_columns(df_t_pr.columns.tolist(), p.DICT_COLS_TOKEEP), errors="raise", inplace=True)
	# df_r_pr.rename(columns=rename_columns(df_r_pr.columns.tolist(), p.DICT_COLS_TOKEEP), errors="raise", inplace=True)

	df_bins_copy = deepcopy(bins)
	df_bins_copy.columns = ["Chromosome", "Start", "End", "len"]
	df_bins_pr = pr.PyRanges(df_bins_copy)

	return df_t_pr, df_r_pr, df_bins_pr


def read_pyranges_new(df_t, df_r, bins):
	"""
	Load data as pyranges
	"""
	df_t_pr = pr.PyRanges(
		df_t.rename(columns=rename_columns(df_t.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))
	df_r_pr = pr.PyRanges(
		df_r.rename(columns=rename_columns(df_r.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))
	df_bins_pr = pr.PyRanges(
		bins.rename(columns=rename_columns(bins.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))

	return df_t_pr, df_r_pr, df_bins_pr


def get_bin_len(s_bins):
	"""
	Get length of each bin
	"""
	return [(s_bins[i], s_bins[i + 1] - s_bins[i]) for i in range(0, len(s_bins) - 1)]


def transform_fragments2breakpoints(df_cycle):
	"""
	Creates the list of breakpoints

	Arguments:
		df_cycle:
	Returns:
		list of tuples
	"""
	breakpoints = pd.DataFrame(columns=["chr1", "start", "chr2", "end", "idx1", "idx2", "strand", "circ_id"])
	circ_id = df_cycle[ht.CIRC_ID].drop_duplicates().tolist()
	circ_len = df_cycle[[ht.CIRC_ID, ht.START]].groupby(ht.CIRC_ID).count().reset_index()
	circ_len.columns = [ht.CIRC_ID, "count"]

	for c in circ_id:
		df_temp = df_cycle[df_cycle[ht.CIRC_ID] == c]
		df_idex = df_temp.index.tolist()
		l = circ_len[circ_len[ht.CIRC_ID] == c].iloc[0, 1]
		for i in range(0, l):
			row1 = df_temp.iloc[i, :]
			row2 = df_temp.iloc[(i + 1) % l, :]

			start = None
			end = None
			strand = None
			chr1 = None
			chr2 = None
			idx1 = None
			idx2 = None

			if row1[ht.STRAND] == "+" and row2[ht.STRAND] == "+":

				if row1[ht.END] < row2[ht.START]:
					chr1 = row1[ht.CHR]
					chr2 = row2[ht.CHR]
					start = row1[ht.END]
					end = row2[ht.START]
					strand = "".join([row1[ht.STRAND], row2[ht.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[ht.CHR]
					chr2 = row1[ht.CHR]
					start = row2[ht.START]
					end = row1[ht.END]
					strand = "".join([row2[ht.STRAND], row1[ht.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]


			elif row1[ht.STRAND] == "+" and row2[ht.STRAND] == "-":

				if row1[ht.END] < row2[ht.END]:
					chr1 = row1[ht.CHR]
					chr2 = row2[ht.CHR]
					start = row1[ht.END]
					end = row2[ht.END]
					strand = "".join([row1[ht.STRAND], row2[ht.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[ht.CHR]
					chr2 = row1[ht.CHR]
					start = row2[ht.END]
					end = row1[ht.END]
					strand = "".join([row2[ht.STRAND], row1[ht.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]

			elif row1[ht.STRAND] == "-" and row2[ht.STRAND] == "-":

				if row1[ht.START] < row2[ht.END]:
					chr1 = row1[ht.CHR]
					chr2 = row2[ht.CHR]
					start = row1[ht.START]
					end = row2[ht.END]
					strand = "".join([row1[ht.STRAND], row2[ht.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[ht.CHR]
					chr2 = row1[ht.CHR]
					start = row2[ht.END]
					end = row1[ht.START]
					strand = "".join([row2[ht.STRAND], row1[ht.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]

			elif row1[ht.STRAND] == "-" and row2[ht.STRAND] == "+":

				if row1[ht.START] < row2[ht.START]:
					chr1 = row1[ht.CHR]
					chr2 = row2[ht.CHR]
					start = row1[ht.START]
					end = row2[ht.START]
					strand = "".join([row1[ht.STRAND], row2[ht.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[ht.CHR]
					chr2 = row1[ht.CHR]
					start = row2[ht.START]
					end = row1[ht.START]
					strand = "".join([row2[ht.STRAND], row1[ht.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]

			breakpoints = breakpoints.append({
				"chr1": chr1,
				"start": start,
				"idx1": idx1,
				"chr2": chr2,
				"end": end,
				"idx2": idx2,
				"strand": strand,
				"circ_id": c}, ignore_index=True)

	return breakpoints


def get_breakpoints_pairs(df_t, df_r):
	"""
	Get the breakpoints for both true and reconstruction

	Arguments:
		df_t (pd.DataFrame):
		df_r (pd.DataFrame):

	Returns:

	"""
	br_t = transform_fragments2breakpoints(df_t)
	br_r = transform_fragments2breakpoints(df_r)

	return br_t, br_r


def get_feature_cn(df_cycles, bins):
	"""
	Get copy-number binned using defined intervals
	"""

	a = BedTool.from_dataframe(bins[[ht.CHR, ht.START, ht.END]])
	df_new = df_cycles[[ht.CHR, ht.START, ht.END, ht.CN]].sort_values(by=[ht.CHR, ht.START, ht.END]).groupby(
		[ht.CHR, ht.START, ht.END]).sum().reset_index()
	b = BedTool.from_dataframe(df_new[[ht.CHR, ht.START, ht.END, ht.CN]])
	o1 = a.intersect(b, wo=True, loj=True).to_dataframe().iloc[:, [0, 1, 2, 6]]
	o1.columns = [ht.CHR, ht.START, ht.END, ht.CN]
	o1.loc[o1[ht.CN] == ".", ht.CN] = 0
	o1[ht.CN] = pd.to_numeric(o1[ht.CN])
	o1 = o1.groupby([ht.CHR, ht.START, ht.END]).sum().reset_index()

	return o1


#
# def get_feature_cn(cycle_fragments, bins):
# 	"""
# 	Get the structure copy-number binned using defined intervals
#
# 	Args:
# 		cycle_fragments (pd.DataFrame): Contains the fragments intervals of the reconstruction
# 										Example
# 										+--------------+-----------+-----------+-----------+-----------+--------------+
# 										| Chromosome   |     Start |       End |   Circ_id |    CN     |  Strand      |
# 										| (category)   |   (int64) |   (int64) |   (int64) |   (int64) | (category)   |
# 										|--------------+-----------+-----------+-----------+-----------+--------------|
# 										| chr1         |      1000 |      2000 |         1 |        10 | +            |
# 										| chr1         |      7000 |      8000 |         1 |        10 | +            |
# 										| chr1         |      1500 |      2100 |         2 |        20 | +            |
# 										| chr1         |       600 |      2500 |         1 |        10 | -            |
# 										| chr1         |      9000 |     11100 |         2 |        20 | -            |
# 										| chr2         |       100 |      2000 |         2 |        20 | +            |
# 										+--------------+-----------+-----------+-----------+-----------+--------------+
# 		bins (pd.DataFrame): Intervals how to bin the genome
# 	"""
# 	START_A = 2
# 	END_A = 3
# 	START_B = 6
# 	END_B = 7
#
# 	# specify which column to keep
# 	cols = cycle_fragments.columns.tolist()
# 	tokeep = []
# 	for c in p.COLS_TOKEEP_SORTED:
# 		found = 0
# 		for c_map in cols:
# 			if c_map in p.DICT_COLS_TOKEEP and p.DICT_COLS_TOKEEP[c_map] == c:
# 				tokeep.append(c_map)
# 				found = 1
# 				break
# 		if found == 0:
# 			raise ValueError("Column name could not be mapped. ",c)
# 	print("columns to keep", tokeep)
#
# 	a = BedTool.from_dataframe(cycle_fragments[tokeep])
# 	b = BedTool.from_dataframe(bins)
#
# 	overlap = b.window(a, w=10).overlap(cols=[START_A, END_A, START_B, END_B])
# 	overlap.columns = ["#chr", "start", "end", "len", "#chr_frag", "start_frag", "stop_frag",
# 					   "circ_id", "cn", "strand","overlap"]
#
# 	overlap = overlap.to_dataframe()
# 	# ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes']
# 	overlap.columns = ["#chr", "start", "end", "len", "#chr_frag", "start_frag", "stop_frag", "circ_id", "cn", "strand", "overlap"]
#
# 	# set to 0 negative overlapping blocks
# 	overlap.loc[overlap["overlap"] < 0, "overlap"] = 0
#
# 	# convert to numeric
# 	cols = ["start", "end", "cn", "overlap"]
# 	overlap[cols] = overlap[cols].apply(pd.to_numeric, errors='coerce')
#
# 	# compute the cummulative coverage
# 	overlap["hit"] = overlap.apply(lambda x: 1 if abs(x['end'] - x['start']) == x["overlap"] else 0, axis=1)
# 	overlap["coverage_norm"] = overlap.apply(lambda x: x['cn'] * x["overlap"], axis=1)
# 	overlap["coverage_mean"] = overlap.apply(lambda x: x['cn'] * x["hit"], axis=1)
# 	return overlap[["#chr", "start", "end", "len", "coverage_norm", "coverage_mean"]].groupby(
# 		by=["#chr", "start", "end", "len"]).agg({'coverage_norm': 'sum', 'coverage_mean': 'sum'}).reset_index()


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


def sigmoid_unmatched(x, a=p.UNMATCHED):
	z = 1 / (1 + np.exp(0.1 * (a - x)))

	# breakpoints unmatch
	if z < 1:
		return True
	# breakpoints match
	return False


def is_unmatched(row1, row2):
	"""
	Check if breakpoints are unmatched (True).
	"""
	if row1[ht.CHR1] != row2[ht.CHR1]:
		return True

	if row1[ht.CHR2] != row2[ht.CHR2]:
		return True

	d = abs(row1[ht.START] - row2[ht.START]) + abs(row1[ht.END] - row2[ht.END])
	if sigmoid_unmatched(d):
		return True

	return False


def create_cost_matrix(br_t, br_r, unmatched, dist=ddt.EUCLIDIAN):
	"""
	Create the cost matrix for the breakpoints pair.
	Set breakpoints with x-x' or y-y' > unmatched as infinity distance

	Arguments:
		br_t (pd.DataFrame): Breakpoints true
		br_r (pd.DataFrame): Breakpoints reconstructed
			Example
					chr1	start	chr2	end	idx1	idx2	strand	circ_id
				0	chr1	5000	chr1	6100	0	1	+-	1
				1	chr1	100		chr1	1000	1	0	-+	1
				2	chr1	6100	chr1	15000	2	3	++	2
				3	chr1	3000	chr1	26000	2	3	++	2
		dist:

	Return:

	"""
	m = np.zeros(shape=(br_t.shape[0], br_r.shape[0]))

	for i in range(0, br_t.shape[0]):
		for j in range(0, br_r.shape[0]):
			m[i, j] = dict_distance_paired_breakpoints[dist](br_t.loc[i, ht.START], br_t.loc[i, ht.END],
															 br_r.loc[j, ht.START], br_r.loc[j, ht.END])

		# # compute distance between breakpoint only for those who are not too far away
		# if not is_unmatched(br_t.loc[i, :], br_r.loc[j, :], unmatched):
		# 	m[i, j] = dict_distance_paired_breakpoints[dist](br_t.loc[i, ht.START], br_t.loc[i, ht.END],
		# 													 br_r.loc[j, ht.START], br_r.loc[j, ht.END])
		# else:
		# 	m[i, j] = p.INF

	print(m)
	return m


def create_bipartite_graph(br_t, br_r, unmatched, dist=ddt.EUCLIDIAN):
	"""
	Create bipartite graph using the breakpoints pairs as nodes and
	edges weighted by the distance (cost) between the pairs).
	Set as unmatched those breakpoints for which x-x' or y-y' > unmatched
	"""
	m = create_cost_matrix(br_t, br_r, unmatched, dist=dist)

	visited_nodes = collections.defaultdict(list)

	true_nodes = {"t" + str(v): v for v in range(0, m.shape[0])}
	reconstruct_nodes = {"r" + str(v): v for v in range(0, m.shape[1])}

	threshold = dict_distance_paired_breakpoints_thresholds[dist]
	list_of_tuples = []

	for k in true_nodes:
		for l in reconstruct_nodes:
			# exclude if distance larger then this
			if m[true_nodes[k], reconstruct_nodes[l]] <= threshold:
				# print()
				# print("included")
				# print(true_nodes[k],reconstruct_nodes[l])
				list_of_tuples.append((k, l, m[true_nodes[k], reconstruct_nodes[l]]))
				visited_nodes[k].append(l)
	# else:
	# 	print()
	# 	print("not included")
	# 	print(true_nodes[k], reconstruct_nodes[l])

	# add minimal amount of edges which will make the graph connected
	max_degree = 0
	max_node = -1
	for k in visited_nodes:
		if len(visited_nodes[k]) > max_degree:
			max_degree = len(visited_nodes[k])
			max_node = k

	# append edges from max node (from true set) to all reconstruct nodes
	for l in reconstruct_nodes:
		if l not in visited_nodes[max_node]:
			list_of_tuples.append((max_node, l, p.MAX_COST))

	G = nx.Graph()
	G.add_weighted_edges_from(list_of_tuples)
	return G, true_nodes, reconstruct_nodes


def find_matching_breakpoints(G, br_t, br_r, t_nodes, r_nodes):
	"""
	Find best matching breakpoints.

	Returns:
		breakpoint_match=[(2,2,"+"),(3,3,"-")...]
	"""
	matches = nx.algorithms.bipartite.minimum_weight_full_matching(G)

	breakpoint_match = []
	# transform this match back to a breakpoint point match
	for k in matches:
		# keep just matches from true to reconstruct
		# (this is symmetrical the other way around)
		if k in t_nodes:
			id_pair_t = t_nodes[k]
			id_pair_r = r_nodes[matches[k]]
			# x1,y1 = br_t.loc[id_pair_t, h.START], br_t.loc[id_pair_t, h.END]
			# x2,y2 = br_r.loc[id_pair_r, h.START], br_r.loc[id_pair_r, h.END]
			breakpoint_match.append((id_pair_t, id_pair_r, G.get_edge_data(k, matches[k])['weight']))
	# breakpoint_match.append((y1, y2))
	return matches, breakpoint_match


def compute_breakpoint_distance(df_t, df_r, unmatched=10000, distance=ddt.EUCLIDIAN):
	"""
	Compute breakpoint-pairs similarity.
	Do not consider breakpoints for which (x-x') > unmatched or (y-y') > unmatched
	Use `distance` to evaluate if the breakpoints are matched or not.

	"""

	# convert fragments to breakpoints
	br_t, br_r = get_breakpoints_pairs(df_t, df_r)

	# create bipartite graph
	G, t_nodes, r_nodes = create_bipartite_graph(br_t,
												 br_r,
												 unmatched=unmatched,
												 dist=distance)

	# get matches
	matches, breakpoint_match = find_matching_breakpoints(G, br_t, br_r, t_nodes, r_nodes)

	# compute jaccard index
	jd = len(breakpoint_match) / (len(br_t) + len(br_r) - len(breakpoint_match))

	return jd, matches, breakpoint_match


def remove_ref2methods(dict_configs):
	"""
	Remove from dictionary all values which are pointers to methods
	"""
	for key in dict_configs[ht.CONFIGS]:
		if key in [ht.BREAKPOINT_DISTANCE, ht.GENOMIC_FOOTPRINT]:
			for d in dict_configs[ht.CONFIGS][key]:
				del dict_configs[ht.CONFIGS][key][d][ht.DEFINITION]
	return dict_configs


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

	# 1. Genome binning based on the breakpoints union

	# as data.frame objects
	df_t, df_r = read_input(t_file, r_file)
	bins = bin_genome(df_t, df_r, margin_size=0)

	# as pyranges objects
	df_t_pr, df_r_pr, df_bins_pr = read_pyranges_new(df_t, df_r, bins)

	# 2. Compute hamming distance and others
	h = get_hamming_score(df_t_pr, df_r_pr, df_bins_pr)
	h_norm = get_hamming_score_norm(df_t_pr, df_r_pr, df_bins_pr)

	dict_metrics[ddt.HAMMING] = h
	dict_metrics[ddt.HAMMING_NORM] = h_norm

	# 3. Compute copy-number similarity
	cv_profile_t = get_feature_cn(df_t, bins)
	cv_profile_r = get_feature_cn(df_r, bins)

	cv_similarity = get_cosine_distance_cn(cv_profile_t, cv_profile_r)
	dict_metrics[ddt.COSINE_DISTANCE] = cv_similarity

	# 4. Breakpoint matching
	jc, matches, breakpoint_matches = compute_breakpoint_distance(df_t, df_r, distance=ddt.RELATIVE_METRIC)
	dict_metrics[ddt.BREAKPOINT_DISTANCE] = jc

	# 5. Penalize for cycles and fragments multiplicity (if the tool decompose in one or more cycles)
	overlap_fragments = get_overlap_fragments_weighted(df_t_pr, df_r_pr, df_bins_pr)
	overlap_cycles = get_overlap_cycles_weighted(df_t_pr, df_r_pr, df_bins_pr)
	dict_metrics[ddt.FRAGMENTS_DISTANCE] = overlap_fragments
	dict_metrics[ddt.CYCLES_DISTANCE] = overlap_cycles

	# 6. Stoichiometry (compute distance of transforming one permutation in the other)

	# 7. Merge results

	# 8. Output
	with open(os.path.join(outdir, 'metrics.json'), 'w', encoding='utf-8') as f:
		final_dict = remove_ref2methods(dict_metrics)
		pprint(final_dict)
		json.dump(final_dict, f, ensure_ascii=False, indent=4, cls=NpEncoder)
	cv_profile_r.to_csv(os.path.join(outdir, 'e2_coverage_profile.txt'), header=True, index=False, sep="\t")
	cv_profile_t.to_csv(os.path.join(outdir, 'e1_coverage_profile.txt'), header=True, index=False, sep="\t")


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
