"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 11:40 AM 08/16/23

Extract genomic features for the cycles set.
"""
import pandas as pd
import numpy as np
import pyranges as pr
from pybedtools import BedTool

from AmpliconComparison.utils.utils import HEADER as ht
from AmpliconComparison.utils.utils import PROPS as p

import warnings
warnings.filterwarnings("ignore")


def get_chromosome_offset(df_t, df_r):
	"""
	Get chromosome offsets for visualization purpose.

	Arguments:
		df_t (pd.DataFrame):
		df_r (pd.DataFrame:

	Returns:
	"""

	df_bins = pd.DataFrame(np.concatenate((df_t[[ht.CHR, ht.START]].values,
										   df_t[[ht.CHR, ht.END]].values,
										   df_r[[ht.CHR, ht.START]].values,
										   df_r[[ht.CHR, ht.END]].values),
										  axis=0))

	df_bins.columns = [ht.CHR, ht.START]
	d3 = df_bins.groupby([ht.CHR]).agg({ht.START: [np.min, np.max]}).reset_index()

	chridx = d3[ht.CHR].drop_duplicates().tolist()

	offsets = {chridx[0]: 0}
	if len(chridx) == 1:
		return offsets

	l = d3.shape[0]
	cumm_offset = 0
	for i in range(1, l):
		max_pos = d3[d3[ht.CHR] == chridx[i - 1]].iloc[0, 1]
		cumm_offset += max_pos
		offsets[chridx[i]] = cumm_offset

	return offsets


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
			 'orientation': 'str',
			 'iscyclic': 'bool'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	# rename columns and check if all columns available
	t_collection.rename(columns=rename_columns(t_collection.columns.tolist(), ht.DICT_HEADER), errors="raise",
						inplace=True)
	r_collection.rename(columns=rename_columns(r_collection.columns.tolist(), ht.DICT_HEADER), errors="raise",
						inplace=True)

	# reorder columns
	keep1 = []
	keep2 = []
	# first check if the columns are in HEADER_SORTED
	for c in ht.HEADER_SORTED:
		if c in t_collection.columns.tolist():
			keep1.append(c)
		if c in r_collection.columns.tolist():
			keep2.append(c)

	for h in t_collection.columns.tolist():
		if h not in ht.HEADER_SORTED:
			keep1.append(h)
	for h in r_collection.columns.tolist():
		if h not in ht.HEADER_SORTED:
			keep2.append(h)

	return t_collection[keep1], r_collection[keep2]


def bin_genome(t_collection, r_collection, margin_size=10000):
	"""
	Bin the intervals by the breakpoints union.
	Warning: Setting a margin size can effect the output of different distances
	"""

	df_bins = pd.DataFrame(np.concatenate((t_collection[[ht.CHR, ht.START]].values,
										   t_collection[[ht.CHR, ht.END]].values,
										   r_collection[[ht.CHR, ht.START]].values,
										   r_collection[[ht.CHR, ht.END]].values),
										  axis=0)).drop_duplicates()

	df_bins.columns = [ht.CHR, ht.START]
	df_bins = df_bins.sort_values(by=[ht.CHR, ht.START])

	# rotate with 1 up the start column
	df_bins_suffix = df_bins.tail(-1)
	df_bins_suffix = pd.concat([df_bins_suffix,df_bins.head(1)], axis=0, ignore_index=True)

	df_bins.reset_index(drop=True, inplace=True)
	df_bins_suffix.reset_index(drop=True, inplace=True)
	df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
	df_bins.columns = [ht.CHR, ht.START, ht.CHRH2, ht.END]
	# keep only rows with same chr and non-negative distance
	df_bins = df_bins[(df_bins[ht.CHR] == df_bins[ht.CHRH2]) & (abs(df_bins[ht.START] - df_bins[ht.END]) >= 0)]

	# # merge intervals
	# df_bins[ht.END] = df_bins.apply(lambda x: x[ht.END] if abs(x[ht.START] - x[ht.END]) < 50 else x[ht.END]-1, axis=1)
	# df_bins = df_bins[(df_bins[ht.END] - df_bins[ht.START]) > 0]
	# a = BedTool.from_dataframe(df_bins[[ht.CHR, ht.START, ht.END]]).sort()
	# df_bins_merged = a.merge(d=0).to_dataframe()
	# df_bins_merged.columns = [ht.CHR, ht.START, ht.END]
	#
	# # merge very small intervals
	# df_bins_merged[ht.LEN] = df_bins_merged[ht.END] - df_bins_merged[ht.START]
	# chrlist = df_bins_merged[ht.CHR].drop_duplicates().tolist()

	df_bins = df_bins[(df_bins[ht.END] - df_bins[ht.START]) > 0]
	df_bins[ht.LEN] = abs(df_bins[ht.END] - df_bins[ht.START])
	chrlist = df_bins[ht.CHR].drop_duplicates().tolist()

	return df_bins[[ht.CHR, ht.START, ht.END, ht.LEN]], chrlist



def read_pyranges(df_t, df_r, bins):
	"""
	Load data as pyranges.

	Args:
		df_t (pd.DataFrame): First dataframe containing reconstructed amplicons
		df_r (pd.DataFrame): Second dataframe containing reconstructed amplicons
		bins (pd.DataFrame): Genomic bins
	"""
	df_t_pr = pr.PyRanges(
		df_t.rename(columns=rename_columns(df_t.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))
	df_r_pr = pr.PyRanges(
		df_r.rename(columns=rename_columns(df_r.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))
	df_bins_pr = pr.PyRanges(
		bins.rename(columns=rename_columns(bins.columns.tolist(), p.DICT_COLS), errors="raise", inplace=False))

	return df_t_pr, df_r_pr, df_bins_pr

def overlap_vector(e1, e2, bins):
	"""
	Overlap bins and the two reconstructed amplicons files

	Args:
		e1 (pr.PyRanges):
		e2 (pr.PyRanges):
		bins (pr.PyRanges):

	Returns:

		A dataframe containing the overlap between the three intervals.
		Columns 'bins,e1,e2' show the overlap count between the sample and Chromosome:Start-End
		Columns 'h1,h2' represent the indicator vectors of 'e1,e2'
		Column 'bin_enabled' is set to 0 if 'e1' or 'e2' does not overlap with the bin

		header: Chromosome	Start	End	bins	e1	e2	h1	h2	bin_enabled

	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()

	# everything > 1 set to 1
	overlaps['h1'] = overlaps['e1'].apply(lambda x: 1 if x >= 1 else 0)
	overlaps['h2'] = overlaps['e2'].apply(lambda x: 1 if x >= 1 else 0)

	# disable bins which are not covered by amplicons
	overlaps[ht.BIN_ENABLED] = 1
	overlaps.loc[(overlaps.h1 == 0) & (overlaps.h2 == 0), ht.BIN_ENABLED] = 0

	return overlaps

def get_bin_len(s_bins):
	"""
	Get length of each bin
	"""
	return [(s_bins[i], s_bins[i + 1] - s_bins[i]) for i in range(0, len(s_bins) - 1)]

def get_feature_cn(df_cycles, bins):
	"""
	Get copy-number binned using defined intervals
	"""

	# bins
	a = BedTool.from_dataframe(bins[[ht.CHR,
									 ht.START,
									 ht.END,
									 ht.BIN_ENABLED]])

	# sum all cn for same bin
	df_new = df_cycles[[ht.CHR, ht.START, ht.END, ht.CN]].\
		groupby([ht.CHR, ht.START, ht.END], observed=True, as_index=False).sum()
	df_new.columns = [ht.CHR, ht.START, ht.END, ht.CN]

	b = BedTool.from_dataframe(df_new[[ht.CHR, ht.START, ht.END, ht.CN]])

	# intersect them
	o1 = a.intersect(b, wo=True, loj=True).to_dataframe().iloc[:, [0, 1, 2, 3, 7]]
	o1.columns = [ht.CHR, ht.START, ht.END, ht.BIN_ENABLED, ht.CN]
	o1.loc[o1[ht.CN] == ".", ht.CN] = 0

	o1[ht.CN] = o1.apply(lambda x: round(float(x[ht.CN]),0), axis=1)
	o1 = o1.groupby([ht.CHR, ht.START, ht.END, ht.BIN_ENABLED], observed=True, as_index=False).sum()

	return o1


