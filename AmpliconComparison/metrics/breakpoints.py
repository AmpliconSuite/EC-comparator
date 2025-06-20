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

from AmpliconComparison.utils.utils import PROPS as p
from AmpliconComparison.utils.utils import DDT as ddt
from AmpliconComparison.utils.utils import HEADER as ht
from AmpliconComparison.metrics.distances import *

import numpy as np
import scipy as sp
import scipy.optimize  # call as sp.optimize
from tabulate import tabulate
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
	breakpoints = pd.DataFrame(
		columns=[
			ht.CHR1,
			ht.START,
			ht.CHR2,
			ht.END,
			ht.IDX1,
			ht.IDX2,
			ht.STRAND,
			ht.CIRC_ID,
			ht.ISCYCLIC,
			ht.CN,
		]
	)
	# get all circ_ids
	circ_id = df_cycle[ht.CIRC_ID].drop_duplicates().tolist()
	# count fragments per cycle
	circ_len = (
		df_cycle[[ht.CIRC_ID, ht.START]].groupby(ht.CIRC_ID).count().reset_index()
	)
	circ_len.columns = [ht.CIRC_ID, "count"]
	# get all the columns
	df_col = df_cycle.columns.tolist()

	for c in circ_id:
		# print("Detect breakpoints for circle", c)
		df_temp = df_cycle[df_cycle[ht.CIRC_ID] == c]
		df_idex = df_temp.index.tolist()
		# count fragments
		count_fragments = circ_len[circ_len[ht.CIRC_ID] == c].iloc[0, 1]

		# by default the amplicon is cyclic
		iscyclic = True

		# if we have the information about cyclic / acyclic amplicon
		if ht.ISCYCLIC in df_col:
			iscyclic = df_temp.loc[:, ht.ISCYCLIC].tolist()[0]

		# print("Iscyclic", iscyclic)

		# single fragment path
		if count_fragments == 1:
			if iscyclic is False:
				# there are no breakpoints
				print("Path id", c, "is a one fragment acylic path")
			else:
				# head to tail
				chr1 = chr2 = df_temp.loc[:, ht.CHR].tolist()[0]
				idx1 = idx2 = df_idex[0]
				strand = "+-"
				start = df_temp.loc[:, ht.END].tolist()[0]
				end = df_temp.loc[:, ht.START].tolist()[0]
				cov = df_temp.loc[:, ht.CN].tolist()[0]

				if start > end:
					tmp = start
					start = end
					end = tmp

				breakpoints = pd.concat(
					[
						breakpoints,
						pd.DataFrame(
							{
								ht.CHR1: chr1,
								ht.START: start,
								ht.IDX1: idx1,
								ht.CHR2: chr2,
								ht.END: end
								+ 10,  # make sure there is 1bp difference between start and end
								ht.IDX2: idx2,
								ht.STRAND: strand,
								ht.CIRC_ID: c,
								ht.ISCYCLIC: iscyclic,
								ht.CN: cov,
							},
							index=[0],
						),
					],
					ignore_index=True,
				)

		# multi-fragment path
		else:
			if iscyclic is False:
				print("Path id", c, "is a multi fragment acylic path")

			for i in range(0, count_fragments):
				row1 = df_temp.iloc[i, :]
				row2 = df_temp.iloc[(i + 1) % count_fragments, :]

				# skip breakpoint if the path is acyclic
				if (i + 1) % count_fragments == 0 and iscyclic is False:
					continue

				start = None
				end = None
				strand = None
				chr1 = None
				chr2 = None
				idx1 = None
				idx2 = None
				cov = row1[ht.CN]

				if row1[ht.STRAND] == "+" and row2[ht.STRAND] == "+":
					# head to tail
					strand = "+-"
					if row1[ht.END] < row2[ht.START]:
						chr1 = row1[ht.CHR]
						chr2 = row2[ht.CHR]
						start = row1[ht.END]
						end = row2[ht.START]
						idx1 = df_idex[i]
						idx2 = df_idex[(i + 1) % count_fragments]
					else:
						chr1 = row2[ht.CHR]
						chr2 = row1[ht.CHR]
						start = row2[ht.START]
						end = row1[ht.END]
						idx2 = df_idex[i]
						idx1 = df_idex[(i + 1) % count_fragments]

				elif row1[ht.STRAND] == "+" and row2[ht.STRAND] == "-":
					# head to head
					strand = "++"
					if row1[ht.END] < row2[ht.END]:
						chr1 = row1[ht.CHR]
						chr2 = row2[ht.CHR]
						start = row1[ht.END]
						end = row2[ht.END]

						idx1 = df_idex[i]
						idx2 = df_idex[(i + 1) % count_fragments]
					else:
						chr1 = row2[ht.CHR]
						chr2 = row1[ht.CHR]
						start = row2[ht.END]
						end = row1[ht.END]
						idx2 = df_idex[i]
						idx1 = df_idex[(i + 1) % count_fragments]

				elif row1[ht.STRAND] == "-" and row2[ht.STRAND] == "-":
					# tail to head however this is the same as head to tail in a circular space
					strand = "+-"
					if row1[ht.START] < row2[ht.END]:
						chr1 = row1[ht.CHR]
						chr2 = row2[ht.CHR]
						start = row1[ht.START]
						end = row2[ht.END]
						idx1 = df_idex[i]
						idx2 = df_idex[(i + 1) % count_fragments]
					else:
						chr1 = row2[ht.CHR]
						chr2 = row1[ht.CHR]
						start = row2[ht.END]
						end = row1[ht.START]
						idx2 = df_idex[i]
						idx1 = df_idex[(i + 1) % count_fragments]

				elif row1[ht.STRAND] == "-" and row2[ht.STRAND] == "+":
					# tail to tail
					strand = "--"
					if row1[ht.START] < row2[ht.START]:
						chr1 = row1[ht.CHR]
						chr2 = row2[ht.CHR]
						start = row1[ht.START]
						end = row2[ht.START]
						idx1 = df_idex[i]
						idx2 = df_idex[(i + 1) % count_fragments]
					else:
						chr1 = row2[ht.CHR]
						chr2 = row1[ht.CHR]
						start = row2[ht.START]
						end = row1[ht.START]
						idx2 = df_idex[i]
						idx1 = df_idex[(i + 1) % count_fragments]

				breakpoints = pd.concat(
					[
						breakpoints,
						pd.DataFrame(
							{
								ht.CHR1: chr1,
								ht.START: start,
								ht.IDX1: idx1,
								ht.CHR2: chr2,
								ht.END: end
								+ 10,  # make sure there is 1bp difference between start and end
								ht.IDX2: idx2,
								ht.STRAND: strand,
								ht.CIRC_ID: c,
								ht.ISCYCLIC: iscyclic,
								ht.CN: cov,
							},
							index=[0],
						),
					],
					ignore_index=True,
				)

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

	if x is None:
		return True

	# [0,1]
	z = 1 / (1 + np.exp(0.8 * (a - x)))

	# breakpoints unmatch
	if z == 1:
		return True
	# breakpoints match
	return False


def f3(x, d=1):
	return np.exp(-x / d)


def f2_sim(x, d=3000, a=3):
	return d / a * np.log(1.0 * d / x)


def f2(x, d=3000, a=3):
	return d * np.exp(-1.0 * a * x / d)


def f1(x):
	return 1.0 * x / (x * x)


def hssp_curve_derivative(x, a=0.7, b=-97, c=8400):
	return -a * c * (x ** (-a - 1)) * (1 + np.exp(-1.0 * x / c)) - x ** (-a) * np.exp(
		-1.0 * x / c
	)


def linear_function(x, d):
	return -x + d


def inverse_function(x, d):
	return 1.0 * d / x


# Define the HSSP function
def hssp_curve(x, a=0.7, b=-97, c=8400):
	# return c * x ** -a * (1 + e ** -x/c) + b
	# a = 0.5
	# b = 0
	# c = 100
	return c * (x ** (-a)) * (1 + np.exp(-x / c)) + b


def is_unmatched_unlinear(row1, row2, match_nonlinear):
	"""
	Use a nonlinear function to match breakpoint-pairs.
	"""
	if (
		row1[ht.CHR1] == row2[ht.CHR1]
		and row1[ht.CHR2] != row2[ht.CHR2]
		or row1[ht.CHR1] == row2[ht.CHR2]
		and row1[ht.CHR2] != row2[ht.CHR1]
	):
		return True

	dx = (
		abs(row1[ht.START] - row2[ht.START])
		if row1[ht.CHR1] == row2[ht.CHR1]
		else abs(row1[ht.START] - row2[ht.END])
	)
	dy = (
		abs(row1[ht.END] - row2[ht.END])
		if row1[ht.CHR1] == row2[ht.CHR1]
		else abs(row1[ht.END] - row2[ht.START])
	)

	if dx > dy:
		tmp = dx
		dx = dy
		dy = tmp

	if match_nonlinear == ddt.MATCH_NONLINEAR_F1:
		# match
		if hssp_curve(dx) > dy:
			return False
		else:
			return True
	elif match_nonlinear == ddt.MATCH_NONLINEAR_F2:
		# match
		if f2(dx) > dy:
			return False
		else:
			return True

	return False


def is_unmatched(
	row1,
	row2,
	unmatched_dist=ddt.MANHATTAN,
	unmatched_threshold=ddt.MANHATTAN_THRESHOLD,
	strandness=True,
):
	"""
	Check if breakpoints are matched or unmatched (True).
	"""
	if strandness and row1[ht.STRAND] != row2[ht.STRAND]:
		return True

	# manhattan distance larger than the sigmoid function
	d = None

	if (
		row1[ht.CHR1] == row2[ht.CHR1]
		and row1[ht.CHR2] != row2[ht.CHR2]
		or row1[ht.CHR1] == row2[ht.CHR2]
		and row1[ht.CHR2] != row2[ht.CHR1]
	):
		return True

	if unmatched_dist == ddt.MANHATTAN:
		if row1[ht.CHR1] == row2[ht.CHR1]:
			d = abs(row1[ht.START] - row2[ht.START]) + abs(row1[ht.END] - row2[ht.END])
		elif row1[ht.CHR1] == row2[ht.CHR2]:
			d = abs(row1[ht.START] - row2[ht.END]) + abs(row1[ht.END] - row2[ht.START])
		if sigmoid_unmatched(d, a=unmatched_threshold):
			return True
	elif unmatched_dist == ddt.EUCLIDIAN:
		if row1[ht.CHR1] == row2[ht.CHR1]:
			d = math.sqrt(
				(row1[ht.START] - row2[ht.START]) ** 2
				+ (row1[ht.END] - row2[ht.END]) ** 2
			)
		elif row1[ht.CHR1] == row2[ht.CHR2]:
			d = math.sqrt(
				(row1[ht.START] - row2[ht.END]) ** 2
				+ (row1[ht.END] - row2[ht.START]) ** 2
			)
		if sigmoid_unmatched(d, a=unmatched_threshold):
			return True

	return False


def create_cost_matrix(br_t, br_r, dist=ddt.EUCLIDIAN):
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
	print("Create cost matrix using distance: ", dist)
	m = np.zeros(shape=(br_t.shape[0], br_r.shape[0]))

	if br_t.shape[0] == 0 or br_r.shape[0] == 0:
		return None

	# compute the distance between the breakpoints-pair
	for i in range(0, br_t.shape[0]):
		for j in range(0, br_r.shape[0]):
			m[i, j] = dict_distance_paired_breakpoints[dist](
				br_t.loc[i, ht.CHR1],
				br_t.loc[i, ht.START],
				br_t.loc[i, ht.CHR2],
				br_t.loc[i, ht.END],
				br_t.loc[i, ht.CN],
				br_r.loc[j, ht.CHR1],
				br_r.loc[j, ht.START],
				br_r.loc[j, ht.CHR2],
				br_r.loc[j, ht.END],
				br_r.loc[j, ht.CN],
			)
	return m


def remove_unmatched(
	m,
	br_t,
	br_r,
	dist=ddt.EUCLIDIAN,
	dist_t=ddt.EUCLIDIAN_THRESHOLD,
	unmatched_dist=ddt.MANHATTAN,
	unmatched_threshold=ddt.MANHATTAN_THRESHOLD,
	strandness=True,
	match_nonlinear=None,
):
	"""
	Set breakpoints with x-x' or y-y' > unmatched as infinity distance
	"""
	# remove edges in the cost matrix if not meeting the criteria
	for i in range(0, br_t.shape[0]):
		for j in range(0, br_r.shape[0]):
			if dist == ddt.EUCLIDIAN and is_unmatched(
				br_t.loc[i, :],
				br_r.loc[j, :],
				unmatched_dist=unmatched_dist,
				unmatched_threshold=unmatched_threshold,
				strandness=strandness,
			):
				m[i, j] = np.inf
			if dist == ddt.EUCLIDIAN and m[i, j] >= dist_t:
				m[i, j] = np.inf
			if dist == ddt.RELATIVE_METRIC and m[i, j] >= dist_t:
				m[i, j] = np.inf

			# elif match_nonlinear:
			# 	if is_unmatched_unlinear(br_t.loc[i, :], br_r.loc[j, :],match_nonlinear, dict_distance_unlinear[match_nonlinear]):
			# 		m[i, j] = np.nan

	# replace all nan values with a max value
	# max_ = max(np.nanmax(m), threshold)
	# new_max = max_ + 1
	# m = np.nan_to_num(m, nan=new_max)
	return m


def create_bipartite_graph(m, br_t, br_r, nodes_weight_function=ddt.COVERAGE):
	"""
	Create bipartite graph using the breakpoints pairs as nodes and
	edges weighted by the distance (cost) between the pairs).
	Set as unmatched those breakpoints for which x-x' or y-y' > unmatched
	"""

	# check if all breakpoints are unmatched
	if np.isnan(m).all():
		return None, {}, {}

	visited_nodes_t = collections.defaultdict(list)
	visited_nodes_r = collections.defaultdict(list)

	# rename the nodes
	true_nodes = {
		"t" + str(v): [v, dict_weight_nodes[nodes_weight_function](br_t, v, m)]
		for v in range(0, m.shape[0])
	}
	reconstruct_nodes = {
		"r" + str(v): [v, dict_weight_nodes[nodes_weight_function](br_r, v, m)]
		for v in range(0, m.shape[1])
	}

	# initialize visited nodes
	list_of_tuples = []

	for k in true_nodes:
		visited_nodes_t[k] = []

	for l in reconstruct_nodes:
		visited_nodes_r[l] = []

	for k in true_nodes:
		for l in reconstruct_nodes:
			# include tuple only if the distance is smaller than threshold
			i = true_nodes[k][0]
			j = reconstruct_nodes[l][0]
			list_of_tuples.append((k, l, m[i, j]))
			visited_nodes_t[k].append(l)
			visited_nodes_r[l].append(k)

	# null nodes which connect the unmatched nodes from the other shore
	reconstruct_nodes[p.RNULL] = [-1, np.inf]
	true_nodes[p.TNULL] = [-1, np.inf]

	# connect null and other node if the node does not have any
	for k in visited_nodes_t:
		# unmatched breakpoint
		if k != p.TNULL:
			list_of_tuples.append((k, p.RNULL, np.inf))

	for k in visited_nodes_r:
		# unmatched breakpoint
		if k != p.RNULL:
			list_of_tuples.append((p.TNULL, k, np.inf))

	list_of_tuples.append((p.TNULL, p.RNULL, np.inf))

	G = nx.Graph()

	# Add weighted edges
	G.add_weighted_edges_from(list_of_tuples)

	# Add weights to nodes
	node_weights = {
		key: value[1] for key, value in {**true_nodes, **reconstruct_nodes}.items()
	}
	nx.set_node_attributes(G, node_weights, "weight")

	return G, true_nodes, reconstruct_nodes


def find_matching_breakpoints(G, t_nodes, r_nodes, threshold_max_value):
	"""
	Find best matching breakpoints.

	Returns:
			breakpoint_match=[(2,2,"+"),(3,3,"-")...]
	"""
	#  - minimum_weight_full_matching specifically wants a full matching with minimum weight: in a bipartite graph with bipartition (𝑈,𝑉)
	# , it wants a matching of size min{|𝑈|,|𝑉|}
	# . If there is no such matching at all, it gives an error.
	# - min_weight_matching (which is slower because it uses the blossom algorithm for matchings in general graphs)
	#   finds the maximum matching with minimum weight: it does not care how big a maximum matching is.

	# print(nx.adjacency_matrix(G))
	# matches = nx.algorithms.bipartite.minimum_weight_full_matching(G)

	# adapt network version
	# left, right = nx.bipartite.sets(G, None)
	# U = list(left)
	# V = list(right)
	U = list(t_nodes.keys())
	V = list(r_nodes.keys())

	weights_sparse = nx.algorithms.bipartite.biadjacency_matrix(
		G, row_order=U, column_order=V, weight="weight", format="coo"
	)

	weights = np.array(weights_sparse.toarray())
	weights[np.isinf(weights)] = threshold_max_value + 1
	left_matches = sp.optimize.linear_sum_assignment(weights)

	# keep only left matches (graph is undirected)
	matches = {U[u]: V[v] for u, v in zip(*left_matches)}

	breakpoint_match = []
	todel = [p.RNULL, p.TNULL]
	# transform this match back to a breakpoint point match
	for k in matches:
		# keep just matches from true to reconstruct
		# (this is symmetrical the other way around)
		if k in [p.RNULL, p.TNULL] or matches[k] in [p.RNULL, p.TNULL]:
			todel.append(k)
		else:
			if k in t_nodes and matches[k] in r_nodes:
				id_pair_t = t_nodes[k]
				id_pair_r = r_nodes[matches[k]]
				props = G.get_edge_data(k, matches[k])
				# print("Pair ids: [s1], [s2] ", id_pair_t, id_pair_r)
				# print(id_pair_t, id_pair_r, props)
				if props is not None and props["weight"] < threshold_max_value:
					breakpoint_match.append(
						(
							# id_pair_t[0],  # id true
							# id_pair_r[0],  # id reconstruct
							k,  # id true
							matches[k],  # id reconstruct
							props["weight"],  # edge weight
							id_pair_t[1],  # node weight
							id_pair_r[1],  # node weight
						)
					)
				else:
					todel.append(k)

	for k in todel:
		if k in matches:
			del matches[k]
	return matches, breakpoint_match


def compute_jc_distance_unweighted(breakpoint_match, t_nodes, r_nodes, G):

	print(
		"Breakpoint maching unweighted: Match score: ",
		2 * len(breakpoint_match),
		", Total score: ",
		(len(t_nodes) + len(r_nodes) - 2),
	)  # -2 excludes the "null" nodes
	return 1 - 1.0 * 2 * len(breakpoint_match) / (len(t_nodes) + len(r_nodes) - 2)


def compute_jc_distance_cn_weighted(breakpoint_match, t_nodes, r_nodes, G):
	match_score = 0
	total_score = 0

	for (
		node1_id,
		node2_id,
		edge_match_weight,
		node1_weight,
		node2_weight,
	) in breakpoint_match:
		match_score += node1_weight + node2_weight

	# weight of the individual
	for key, value in {**t_nodes, **r_nodes}.items():

		# key - node id
		# value - tuple (breakpoint id, weight)
		if key not in {p.RNULL, p.TNULL}:
			total_score += value[1]

	print(
		"Breakpoint maching cn weighted: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)
	return max(0, 1 - 1.0 * match_score / total_score)


def avg_connections(node, G):

	# print("Neigbors for: ", node)
	# print(G.neighbors(node))
	neighbors_weights = [
		G.nodes[nbr].get("weight", None)
		for nbr in G.neighbors(node)
		if G[node][nbr]["weight"] < np.inf
	]
	# print(neighbors_weights)
	return np.mean(neighbors_weights)


def compute_jc_distance_weighted_avg(breakpoint_match, t_nodes, r_nodes, G):

	match_score = 0
	total_score = 0

	# print()
	# print(breakpoint_match)
	# m = nx.to_numpy_array(G, weight="weight")
	# print(tabulate(m, tablefmt="grid", floatfmt=".3f"))

	for (
		node1_id,
		node2_id,
		edge_match_weight,
		node1_weight,
		node2_weight,
	) in breakpoint_match:
		wn1 = avg_connections(node2_id, G)
		wn2 = avg_connections(node1_id, G)
		match_score += wn1 + wn2

	# weight of the individual
	for key, value in {**t_nodes, **r_nodes}.items():

		# key - node id
		# value - tuple (breakpoint id, weight)
		if key not in {p.RNULL, p.TNULL}:
			total_score += value[1]

	print(
		"Breakpoint maching cn avg weighted: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)
	return 1 - 1.0 * match_score / total_score


def compute_jc_distance_unweighted_confidence(breakpoint_match, t_nodes, r_nodes, G):
	"""
	Compute the breakpoint distance by including the gaussian confidence per breakpoint match.
	JD = 1 - 2*(0.529+0.9615)/(1+1+1+1+1+1)= 0.50, with 3 bp for s1 and 3 bp for s2 and 2 breakpoint-pairs matching.
	"""
	# b[2] is the edge weight and corresponds to the gaussian distance between the breakpoints
	# 1 - b[2] is the gaussian confidence
	match_score = 2 * np.sum(np.array([1 - b[2] for b in breakpoint_match]))
	total_score = len(t_nodes) + len(r_nodes) - 2
	print(
		"Breakpoint maching unweighted + gaussian confidence: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)  # -2 excludes the "null" nodes
	return 1 - 1.0 * match_score / total_score


def compute_jc_distance_cn_weighted_confidence(breakpoint_match, t_nodes, r_nodes, G):
	"""
	Compute the breakpoint distance weighted by the coverage,
	by including the gaussian confidence per breakpoint match.

	JD = 1 - (10*0.529+100*0.529+100*0.9615+10*0.9615)/(10+100+100+10+90+10)=0.48,
	with 3 bp for s1 and 3 bp for s2 and 2 breakpoint-pairs matching.
	"""
	match_score = 0
	total_score = 0

	for (
		node1_id,
		node2_id,
		edge_match_weight,
		node1_weight,
		node2_weight,
	) in breakpoint_match:
		gaussian_confidence = 1 - edge_match_weight
		match_score += gaussian_confidence * (node1_weight + node2_weight)

	# weight of the individual
	for key, value in {**t_nodes, **r_nodes}.items():
		# key - node id
		# value - tuple (breakpoint id, weight)
		if key not in {p.RNULL, p.TNULL}:
			total_score += value[1]

	print(
		"Breakpoint maching cn weighted  + gaussian confidence: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)
	return max(0, 1 - 1.0 * match_score / total_score)


def compute_jc_distance_cn_weighted_avg_confidence(
	breakpoint_match, t_nodes, r_nodes, G
):
	"""
	Compute the breakpoint distance weighted by the coverage average,
	by including the gaussian confidence per breakpoint match.

	JD = 1 - (50*0.529+100*0.529+100*0.9615+10*0.9615)/(10+100+100+10+90+10)=0.42,
	with 3 bp for s1 and 3 bp for s2 and 2 breakpoint-pairs matching.
	"""
	match_score = 0
	total_score = 0

	for (
		node1_id,
		node2_id,
		edge_match_weight,
		node1_weight,
		node2_weight,
	) in breakpoint_match:
		wn1 = avg_connections(node2_id, G)
		wn2 = avg_connections(node1_id, G)
		gaussian_confidence = 1 - edge_match_weight
		match_score += gaussian_confidence * (wn1 + wn2)

	# weight of the individual
	for key, value in {**t_nodes, **r_nodes}.items():

		# key - node id
		# value - tuple (breakpoint id, weight)
		if key not in {p.RNULL, p.TNULL}:
			total_score += value[1]

	print(
		"Breakpoint maching cn avg weighted + gaussian confidence: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)
	return max(0, 1 - 1.0 * match_score / total_score)


def compute_bp_gaussian_distance(m, br_t, br_r, how):
	"""
	Compute breakpoint distance without breakpoint matching, instead use only the total gaussian contribution.
	JD = 1 - 3,07/(1+1+1+1+1+1)= 0.38 (unweighted)
	JD = 1 - (10*0.530+100*0.8905+100*0.9615+10*0.9625+90*0.3615+10*0)/(10+100+100+10+90+10)= 0.27 (weighted)
	"""
	m_test = deepcopy(m)
	m_test[m_test == np.inf] = 0
	m_test[m_test != 0] = 1 - m_test[m_test != 0]
 
	print(tabulate(m_test, tablefmt="grid", floatfmt=".3f"))

	row_sums = m_test.sum(axis=1)  
	col_sums = m_test.sum(axis=0)
 
	match_score = None
	total_score = None
	
	if how == ddt.BP_GAUSSIAN_CONFIDENCE_UNWEIGHTED:
		match_score = np.sum(row_sums) + np.sum(col_sums)
		total_score = m_test.shape[0] + m_test.shape[1]
	elif how == ddt.BP_GAUSSIAN_CONFIDENCE_CN_WEIGHTED:
		weights = [row[9] for row in br_t.itertuples(index=False)]
		weights += [row[9] for row in br_r.itertuples(index=False)]
	
		rows_weighted = [x * weights[i] for i,x in enumerate(row_sums)]
		cols_weighted = [x * weights[br_t.shape[0] + i] for i,x in enumerate(col_sums)]

		match_score = np.sum(rows_weighted) + np.sum(cols_weighted)
		total_score = np.sum(weights)
	else:
		raise Exception(
			f"{how} not defined! Please use one of the following options for breakpoint distance calculation: {ddt.BP_GAUSSIAN_CONFIDENCE_UNWEIGHTED} or {ddt.BP_GAUSSIAN_CONFIDENCE_CN_WEIGHTED}"
		)
	 
	print(
		f"{how}: Match score: ",
		match_score,
		", Total score: ",
		total_score,
	)
	return 1 - 1.0 * match_score / total_score


def compute_jc_distance(
	breakpoint_match, t_nodes, r_nodes, G, how=ddt.BP_MATCH_UNWEIGHTED
):
	"""
	Compute jaccard distance between the matched and unmached breakpoints

	Arguments:
			breakpoint_match (list): List of tuples (node1_id, node2_id, edge_match_weight, node1_weight, node2_weight)
			t_nodes (dict): Dictionary of all breakpoint pairs in first sample (node_id, node_weight)
			r_nodes (dict): Dictionary of all breakpoint pairs in second sample (node_id, node_weight)
			how (str): How you do the matching

	"""
	print("Compute JD using: ", how)
	jd = dict_breakpoint_distance_calculation[how](
		breakpoint_match, t_nodes, r_nodes, G
	)
	return jd


def compute_breakpoint_distance(
	br_t,
	br_r,
	distance=ddt.EUCLIDIAN,
	distance_threshold=ddt.EUCLIDIAN_THRESHOLD,
	unmatched_dist=ddt.UNMATCHED_DISTANCE,
	unmatched_threshold=ddt.UNMATCHED_THRESHOLD,
	how=ddt.BP_MATCH_UNWEIGHTED,
	match_nonlinear=None,
):
	"""
	Compute breakpoint-pairs distance.
	Do not consider breakpoints for which |x-x'| + |y-y'| >= unmatched
	Use `distance` to evaluate if the breakpoints are matched or not.
	USe `how` to compute the breakpoint-pairs distance.
	"""

	# 1. create cost matrix
	m = create_cost_matrix(br_t, br_r, dist=distance)

	# 2. remove unmatched edges
	if unmatched_dist in [ddt.EUCLIDIAN, ddt.MANHATTAN, ddt.RELATIVE_METRIC]:
		m = remove_unmatched(
			m,
			br_t,
			br_r,
			dist=distance,
			dist_t=distance_threshold,
			unmatched_dist=unmatched_dist,
			unmatched_threshold=unmatched_threshold,
			match_nonlinear=match_nonlinear,
		)
	else:
		# if gaussian skip this step
		print("Skip unmatched step because Gaussian distance set.")

	# 3. create bipartite graph
	G, t_nodes, r_nodes = create_bipartite_graph(
		m, br_t, br_r, nodes_weight_function=ddt.COVERAGE
	)

	if G is None:
		return 1, {}, []

	# 4. get matches
	matches, breakpoint_match = find_matching_breakpoints(
		G, t_nodes, r_nodes, threshold_max_value=distance_threshold
	)
 
	# 5. breakpoint distance
	if how in [
		ddt.BP_GAUSSIAN_CONFIDENCE_UNWEIGHTED,
		ddt.BP_GAUSSIAN_CONFIDENCE_CN_WEIGHTED,
	]:
		# 5. gaussian distance
		jd = compute_bp_gaussian_distance(m, br_t, br_r, how=how)
	else:
		# 5. compute jaccard distance
		jd = compute_jc_distance(breakpoint_match, t_nodes, r_nodes, G, how=how)

	return round(jd, 3), matches, breakpoint_match, G


dict_breakpoint_distance_calculation = {
	ddt.BP_MATCH_UNWEIGHTED: compute_jc_distance_unweighted,
	ddt.BP_MATCH_CN_WEIGHTED: compute_jc_distance_cn_weighted,
	ddt.BP_MATCH_CN_WEIGHTED_AVG: compute_jc_distance_weighted_avg,
	ddt.BP_MATCH_UNWEIGHTED_CONFIDENCE: compute_jc_distance_unweighted_confidence,
	ddt.BP_MATCH_CN_WEIGHTED_CONFIDENCE: compute_jc_distance_cn_weighted_confidence,
	ddt.BP_MATCH_CN_WEIGHTED_AVG_CONFIDENCE: compute_jc_distance_cn_weighted_avg_confidence,
}

dict_distance_unlinear = {
	ddt.MATCH_NONLINEAR_F1: ddt.MATCH_NONLINEAR_F1_THRESHOLD,
	ddt.MATCH_NONLINEAR_F2: ddt.MATCH_NONLINEAR_F2_THRESHOLD,
}
