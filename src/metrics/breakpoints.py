"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 11:50 AM 08/16/23

Breakpoint matching logic.
"""
import collections
import pandas as pd
import networkx as nx
import numpy as np

from src.utils.utils import PROPS as p
from src.utils.utils import DDT as ddt
from src.utils.utils import HEADER as ht
from src.metrics.distances import *

import numpy as np
import scipy as sp
import scipy.optimize  # call as sp.optimize

import warnings
warnings.filterwarnings("ignore")


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


def sigmoid_unmatched(x, a=ddt.SIGMOID_THRESHOLD):
	z = 1 / (1 + np.exp(0.8 * (a - x)))

	# breakpoints unmatch
	if z >= 1:
		return True
	# breakpoints match
	return False


def is_unmatched(row1, row2, strandness=True):
	"""
	Check if breakpoints are unmatched (True).
	"""
	# breakpoints do not belong to same chr
	if row1[ht.CHR1] != row2[ht.CHR1]:
		return True

	if row1[ht.CHR2] != row2[ht.CHR2]:
		return True

	# ensure strandness of the breakpoints
	if strandness and row1[ht.STRAND] != row2[ht.STRAND]:
		return True

	# manhattan distance larger than the sigmoid function
	d = abs(row1[ht.START] - row2[ht.START]) + abs(row1[ht.END] - row2[ht.END])
	if sigmoid_unmatched(d):
		return True

	return False


def create_cost_matrix(br_t, br_r, dist=ddt.EUCLIDIAN, strandness=True):
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

	# compute the distance between the breakpoints-pair
	for i in range(0, br_t.shape[0]):
		for j in range(0, br_r.shape[0]):
			# set distance to infinity for those breakpoints which are far away
			if not is_unmatched(br_t.loc[i, :], br_r.loc[j, :], strandness=strandness):
				m[i, j] = dict_distance_paired_breakpoints[dist](br_t.loc[i, ht.START], br_t.loc[i, ht.END],
																 br_r.loc[j, ht.START], br_r.loc[j, ht.END])
			else:
				m[i, j] = np.nan

	# replace all nan values with a max value
	max = np.nanmax(m)
	new_max = max + 1
	m = np.nan_to_num(m, nan=new_max)
	return m, new_max


def create_bipartite_graph(br_t, br_r, dist, threshold):
	"""
	Create bipartite graph using the breakpoints pairs as nodes and
	edges weighted by the distance (cost) between the pairs).
	Set as unmatched those breakpoints for which x-x' or y-y' > unmatched
	"""
	m, max_value = create_cost_matrix(br_t, br_r, dist=dist)

	visited_nodes_t = collections.defaultdict(list)
	visited_nodes_r = collections.defaultdict(list)

	# rename the nodes
	true_nodes = {"t" + str(v): v for v in range(0, m.shape[0])}
	reconstruct_nodes = {"r" + str(v): v for v in range(0, m.shape[1])}

	# null nodes which connect the unmatched nodes from the other shore
	reconstruct_nodes[p.RNULL] =  -1
	true_nodes[p.TNULL] = -1

	list_of_tuples = []

	for k in true_nodes:
		visited_nodes_t[k] = []

	for l in reconstruct_nodes:
		visited_nodes_r[l] = []

	for k in true_nodes:
		for l in reconstruct_nodes:
			# exclude if distance larger then this
			if m[true_nodes[k], reconstruct_nodes[l]] <= threshold:
				# print()
				# print("included")
				# print(true_nodes[k],reconstruct_nodes[l])
				list_of_tuples.append((k, l, m[true_nodes[k], reconstruct_nodes[l]]))
				visited_nodes_t[k].append(l)
				visited_nodes_r[l].append(k)

	# connect null and other node if the node does not have any
	for k in visited_nodes_t:
		# unmatched breakpoint
		if k != p.TNULL and len(visited_nodes_t[k]) == 0:
			list_of_tuples.append((k, p.RNULL, max_value))

	for k in visited_nodes_r:
		# unmatched breakpoint
		if k != p.RNULL and len(visited_nodes_r[k]) == 0:
			list_of_tuples.append((p.TNULL,k, max_value))

	# add minimal amount of edges which will make the graph connected
	max_degree = 0
	max_node = -1
	for k in visited_nodes_t:
		if k not in [p.RNULL, p.TNULL] and len(visited_nodes_t[k]) > max_degree:
			max_degree = len(visited_nodes_t[k])
			max_node = k

	# append edges from max node (from true set) to all reconstruct nodes
	for l in reconstruct_nodes:
		if l not in visited_nodes_t[max_node]:
			list_of_tuples.append((max_node, l, p.INF))

	G = nx.Graph()
	G.add_weighted_edges_from(list_of_tuples)
	return G, true_nodes, reconstruct_nodes, m

def find_matching_breakpoints(G, t_nodes, r_nodes, cost_matrix):
	"""
	Find best matching breakpoints.

	Returns:
		breakpoint_match=[(2,2,"+"),(3,3,"-")...]
	"""
	#  - minimum_weight_full_matching specifically wants a full matching with minimum weight: in a bipartite graph with bipartition (𝑈,𝑉)
	# , it wants a matching of size min{|𝑈|,|𝑉|}
	# . If there is no such matching at all, it gives an error.
    # - min_weight_matching (which is slower because it uses the blossom algorithm for matchings in general graphs) finds a maximum matching with minimum weight: it does not care how big a maximum matching is.

	matches = nx.algorithms.bipartite.minimum_weight_full_matching(G)

	breakpoint_match = []
	todel = [p.RNULL, p.TNULL]
	# transform this match back to a breakpoint point match
	for k in matches:
		# keep just matches from true to reconstruct
		# (this is symmetrical the other way around)
		if k not in [p.RNULL, p.TNULL] and k in t_nodes and matches[k] not in [p.RNULL, p.TNULL]:
			id_pair_t = t_nodes[k]
			id_pair_r = r_nodes[matches[k]]
			breakpoint_match.append((id_pair_t, id_pair_r,
									 G.get_edge_data(k, matches[k])['weight']))
		else:
			# ignore nodes conneted to RNUL or TNULL
			if matches[k] in [p.RNULL, p.TNULL]:
				todel.append(k)

	for k in todel:
		if k in matches:
			del matches[k]

	return matches, breakpoint_match


def compute_breakpoint_distance(df_t, df_r, distance, threshold):
	"""
	Compute breakpoint-pairs similarity.
	Do not consider breakpoints for which (x-x') > unmatched or (y-y') > unmatched
	Use `distance` to evaluate if the breakpoints are matched or not.

	"""

	# convert fragments to breakpoints
	br_t, br_r = get_breakpoints_pairs(df_t, df_r)

	# create bipartite graph
	G, t_nodes, r_nodes, cost_matrix = create_bipartite_graph(br_t,
												 br_r,
												 dist=distance,
												 threshold=threshold)

	# get matches
	matches, breakpoint_match = find_matching_breakpoints(G, t_nodes, r_nodes, cost_matrix)

	# compute jaccard distance
	jd = 1 - len(breakpoint_match) / (len(br_t) + len(br_r) - len(breakpoint_match))

	return jd, matches, breakpoint_match
