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

	# print(m)
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

	# compute jaccard distance
	jd = 1 - len(breakpoint_match) / (len(br_t) + len(br_r) - len(breakpoint_match))

	return jd, matches, breakpoint_match
