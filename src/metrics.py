import collections
import math
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
from pybedtools import BedTool
import pyranges as pr
import networkx as nx

from utils import HEADER as h
from utils import PROPS as p

def read_input(t_file, r_file):
	"""
	Read input true and reconstruct file
	"""
	dtype = {'#chr': 'str',
			 'start': 'int',
			 'end': 'int',
			 'circ_id': 'str',
			 'estimated_cn': 'double',
			 'strand': 'str'}

	t_collection = pd.read_csv(t_file, header=0, sep="\t", dtype=dtype)
	r_collection = pd.read_csv(r_file, header=0, sep="\t", dtype=dtype)

	return t_collection, r_collection


def bin_genome(t_collection, r_collection, margin_size=10000):
	"""
	Bin the intervals by the breakpoints union.
	"""
	df_bins = pd.DataFrame(np.concatenate((t_collection[["#chr", "start"]].values,
										   t_collection[["#chr", "end"]].values,
										   r_collection[["#chr", "start"]].values,
										   r_collection[["#chr", "end"]].values),
										  axis=0))

	df_bins.columns = ["#chr", "start"]
	df_bins_gr = df_bins.groupby(["#chr"]).agg({"start": [np.min, np.max]}).reset_index()
	df_bins_gr.columns = ["#chr", "start_amin", "start_amax"]
	df_bins_gr["start_amin"] = df_bins_gr.apply(lambda x: max(0, x["start_amin"] - margin_size), axis=1)
	df_bins_gr["start_amax"] = df_bins_gr.apply(lambda x: x["start_amax"] + margin_size, axis=1)

	df_bins = pd.concat([df_bins, df_bins_gr[["#chr", "start_amin"]].rename(columns={"start_amin": 'start'})],
						ignore_index=True)
	df_bins = pd.concat([df_bins, df_bins_gr[["#chr", "start_amax"]].rename(columns={"start_amax": 'start'})],
						ignore_index=True)
	df_bins = df_bins.sort_values(by=["#chr", "start"]).drop_duplicates()

	# rotate with 1 up the start column
	df_bins_suffix = df_bins.tail(-1)
	df_bins_suffix = df_bins_suffix.append((df_bins.head(1)), ignore_index=True)

	df_bins.reset_index(drop=True, inplace=True)
	df_bins_suffix.reset_index(drop=True, inplace=True)
	df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
	df_bins.columns = ["#chr", "start", "#chr2", "end"]
	df_bins["len"] = df_bins["end"].astype(int) - df_bins["start"].astype(int)

	# keep only rows with same chr and positive distance
	df_bins = df_bins[(df_bins["#chr"] == df_bins["#chr2"]) & (df_bins["len"] > 0)]

	return df_bins[["#chr", "start", "end", "len"]]


def get_bin_len(s_bins):
	"""
	Get length of each bin
	"""
	return [(s_bins[i], s_bins[i + 1] - s_bins[i]) for i in range(0, len(s_bins) - 1)]


def get_breakpoints(df_cycle):
	"""
	Creates the list of breakpoints

	Arguments:
		df_cycle:
	Returns:
		list of tuples
	"""
	breakpoints = pd.DataFrame(columns=["chr1", "start", "chr2", "end", "idx1", "idx2", "strand", "circ_id"])
	circ_id = df_cycle[h.CIRC_ID].drop_duplicates().tolist()
	circ_len = df_cycle[[h.CIRC_ID, h.START]].groupby(h.CIRC_ID).count().reset_index()
	circ_len.columns = ["circ_id", "count"]

	for c in circ_id:
		df_temp = df_cycle[df_cycle[h.CIRC_ID] == c]
		df_idex = df_temp.index.tolist()
		l = circ_len[circ_len["circ_id"] == c].iloc[0, 1]
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

			if row1[h.STRAND] == "+" and row2[h.STRAND] == "+":

				if row1[h.END] < row2[h.START]:
					chr1 = row1[h.CHR]
					chr2 = row2[h.CHR]
					start = row1[h.END]
					end = row2[h.START]
					strand = "".join([row1[h.STRAND], row2[h.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[h.CHR]
					chr2 = row1[h.CHR]
					start = row2[h.START]
					end = row1[h.END]
					strand = "".join([row2[h.STRAND], row1[h.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]


			elif row1[h.STRAND] == "+" and row2[h.STRAND] == "-":

				if row1[h.END] < row2[h.END]:
					chr1 = row1[h.CHR]
					chr2 = row2[h.CHR]
					start = row1[h.END]
					end = row2[h.END]
					strand = "".join([row1[h.STRAND], row2[h.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[h.CHR]
					chr2 = row1[h.CHR]
					start = row2[h.END]
					end = row1[h.END]
					strand = "".join([row2[h.STRAND], row1[h.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]

			elif row1[h.STRAND] == "-" and row2[h.STRAND] == "-":

				if row1[h.START] < row2[h.END]:
					chr1 = row1[h.CHR]
					chr2 = row2[h.CHR]
					start = row1[h.START]
					end = row2[h.END]
					strand = "".join([row1[h.STRAND], row2[h.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[h.CHR]
					chr2 = row1[h.CHR]
					start = row2[h.END]
					end = row1[h.START]
					strand = "".join([row2[h.STRAND], row1[h.STRAND]])
					idx2 = df_idex[i]
					idx1 = df_idex[(i + 1) % l]

			elif row1[h.STRAND] == "-" and row2[h.STRAND] == "+":

				if row1[h.START] < row2[h.START]:
					chr1 = row1[h.CHR]
					chr2 = row2[h.CHR]
					start = row1[h.START]
					end = row2[h.START]
					strand = "".join([row1[h.STRAND], row2[h.STRAND]])
					idx1 = df_idex[i]
					idx2 = df_idex[(i + 1) % l]
				else:
					chr1 = row2[h.CHR]
					chr2 = row1[h.CHR]
					start = row2[h.START]
					end = row1[h.START]
					strand = "".join([row2[h.STRAND], row1[h.STRAND]])
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
	br_t = get_breakpoints(df_t)
	br_r = get_breakpoints(df_r)

	return br_t, br_r


def get_feature_cn(cycle_fragments, bins):
	"""
	Get the structure copy-number binned using defined intervals

	Arguments:
		cycle_fragments (pd.DataFrame): Contains the fragments intervals of the reconstruction
										Example
										+--------------+-----------+-----------+-----------+-----------+--------------+
										| Chromosome   |     Start |       End |   Circ_id |    CN     |  Strand      |
										| (category)   |   (int64) |   (int64) |   (int64) |   (int64) | (category)   |
										|--------------+-----------+-----------+-----------+-----------+--------------|
										| chr1         |      1000 |      2000 |         1 |        10 | +            |
										| chr1         |      7000 |      8000 |         1 |        10 | +            |
										| chr1         |      1500 |      2100 |         2 |        20 | +            |
										| chr1         |       600 |      2500 |         1 |        10 | -            |
										| chr1         |      9000 |     11100 |         2 |        20 | -            |
										| chr2         |       100 |      2000 |         2 |        20 | +            |
										+--------------+-----------+-----------+-----------+-----------+--------------+
		bins (pd.DataFrame): Intervals how to bin the genome
	"""
	START_A = 2
	END_A = 3
	START_B = 6
	END_B = 7
	a = BedTool.from_dataframe(cycle_fragments)
	b = BedTool.from_dataframe(bins)

	overlap = b.window(a, w=10).overlap(cols=[START_A, END_A, START_B, END_B])
	overlap.columns = ["#chr", "start", "end", "len", "#chr_frag", "start_frag", "stop_frag", "circ_id", "cn", "strand",
					   "overlap"]

	overlap = overlap.to_dataframe()
	# ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes']
	overlap.columns = ["#chr", "start", "end", "len", "#chr_frag", "start_frag", "stop_frag", "circ_id", "cn", "strand",
					   "overlap"]

	# convert to numeric
	cols = ["start", "end", "cn", "overlap"]
	overlap[cols] = overlap[cols].apply(pd.to_numeric, errors='coerce')

	# compute the cummulative coverage
	overlap["hit"] = overlap.apply(lambda x: 1 if abs(x['end'] - x['start']) == x["overlap"] else 0, axis=1)
	overlap["coverage_norm"] = overlap.apply(lambda x: x['cn'] * x["overlap"], axis=1)
	overlap["coverage_mean"] = overlap.apply(lambda x: x['cn'] * x["hit"], axis=1)
	return overlap[["#chr", "start", "end", "len", "coverage_norm", "coverage_mean"]].groupby(
		by=["#chr", "start", "end", "len"]).agg({'coverage_norm': 'sum', 'coverage_mean': 'sum'}).reset_index()


def get_cos_similarity_cn(e1, e2):
	"""
	Arguments:
		e1 (pd.DataFrame): First reconstruction binned by union of all bins
		e2 (pd.DataFrame): Second reconstruction binned by union of all bins
	"""
	pass

def get_hamming_score(e1, e2, bins):
	"""
	Compute hamming distance between the two reconstructions

	Arguments:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DataFrame): Bins of the intervals union e1 and e2
	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()
	overlaps['hamming'] = overlaps.apply(lambda x: x.e1 ^ x.e2, axis=1)
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
	overlaps['hamming'] = overlaps.apply(lambda x: x.e1 ^ x.e2, axis=1)
	overlaps["len"] = overlaps["End"] - overlaps["Start"]
	overlaps["prod"] = overlaps["hamming"] * overlaps["len"]
	return overlaps["prod"].sum() / overlaps["len"].sum()


def get_overlap_score_weighted(e1, e2, bins):
	"""
	Compute the overlap distance between the two reconstructions

	Args:
		e1 (pd.DataFrame): First reconstruction
		e2 (pd.DataFrame): Second reconstruction
		bins (pd.DatamFrame): Bins of the intervals union e1 and e2
	"""
	grs = {n: s for n, s in zip(["bins", "e1", "e2"], [bins, e1, e2])}
	# Chromosome	Start	End	bins	e1	e2
	overlaps = pr.count_overlaps(grs)
	overlaps = overlaps.as_df()
	overlaps['overlapping_score'] = overlaps.apply(lambda x: abs(x.e1 - x.e2), axis=1)
	overlaps["prod"] = overlaps["overlapping_score"] * (overlaps["End"] - overlaps["Start"])
	return overlaps["prod"].sum()


def euclidian_distance(a, b, x, y):
	return math.sqrt((a - x) ** 2 + (b - y) ** 2)


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

	v1 = np.array([a - x,
				   a - b,
				   a - y])
	v2 = np.array([x - a,
				   x - b,
				   x - y])
	v3 = np.array([b - a,
				   b - x,
				   b - y])
	v4 = np.array([y - a,
				   y - x,
				   y - b])
	### the implementation of scipy.spatial.distance.cosine is
	### 1 - dot(A,B)/||A||*||B|| so this is something we dont want
	### we revert results by subtracting 1
	cos1 =  (1-cosine(v1, v2))
	cos2 = (1-cosine(v3, v4))

	# return (cos1 + cos2) / 2
	return 1 - abs((cos1 + cos2) / 2)


def create_cost_matrix(br_t, br_r, dist="euclidian"):
	"""
	Create the cost matrix for the breakpoints pair.

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
			m[i, j] = dict_distance_paired_breakpoints[dist](br_t.loc[i, h.START], br_t.loc[i, h.END],
															 br_r.loc[j, h.START], br_r.loc[j, h.END])

	return m


def create_bipartite_graph(br_t, br_r, dist="euclidian"):
	"""
	Create bipartite graph using the breakpoints pairs as nodes and
	edges weighted by the distance (cost) between the pairs)
	"""
	m = create_cost_matrix(br_t, br_r, dist=dist)
	visited_nodes = collections.defaultdict(list)

	true_nodes = {"t"+str(v):v for v in range(0, m.shape[0])}
	reconstruct_nodes = {"r" + str(v): v for v in range(0, m.shape[1])}

	threshold = dict_distance_paired_breakpoints_thresholds[dist]
	list_of_tuples = []

	for k in true_nodes:
		for l in reconstruct_nodes:
			# exclude if distance larger then this
			if m[true_nodes[k],reconstruct_nodes[l]] <= threshold:

				# print()
				# print("included")
				# print(true_nodes[k],reconstruct_nodes[l])
				list_of_tuples.append((k,l,m[true_nodes[k],reconstruct_nodes[l]]))
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
			breakpoint_match.append((id_pair_t,id_pair_r, G.get_edge_data(k, matches[k])['weight']))
			# breakpoint_match.append((y1, y2))
	return matches, breakpoint_match




def d1(t_collection, r_collection, sorted_bins, hamming=False):
	pass

def compare_cycles(t_file, r_file, hamming=True, viz=False):
	"""
	Entrypoint: compare the distance between two cycle sets.
	Args:
		t_file (str): First cycle set
		r_file (str): Second cycle set
		hamming (bool): Include copy-number similarity for the total distance if set to True or Hamming distance is copy-number not available
		viz (bool): Enable visualization of the comparison

	Returns:
		Dictionary with different distances
	"""

	# 1. Genome binning based on the breakpoints union

	# 2. Compute hamming distance and others

	# 3. Compute copy-number

	# 4. Breakpoint matching

	# 5. Penalize for cycles multiplicity (if the tool decompose in one or more cycles)

	# 6. Stoichiometry (compute distance of transforming one permutation in the other)

	# 7. Merge results




# distance for the copy-number profile
dict_distance_function = {"hamming": get_hamming_score,
						  "hamming_norm": get_hamming_score_norm,
						  "overlap": get_overlap_score_weighted}

# distance between paired breapoints
dict_distance_paired_breakpoints = {
	"euclidian": euclidian_distance,
	# "auc_triangle": auc_triangle,
	# "auc_trapeze": auc_trapeze,
	# "match_angle": match_angle,
	"match_score": match_score
}

dict_distance_paired_breakpoints_thresholds = {
	"euclidian": 3000,
	"auc_triangle": 500000,
	"auc_trapeze": 50000,
	"match_angle": 1.2,
	"match_score": 0.3
}
