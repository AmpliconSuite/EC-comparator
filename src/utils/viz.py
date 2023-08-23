import math
import pandas as pd
import numpy as np
from numpy import sin, cos, pi, linspace

import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt

# custom module
from src.utils.utils import DDT
from src.utils.utils import HEADER as ht
from src.utils.utils import PROPS

tokeep = [DDT.HAMMING_NORM, DDT.COSINE_DISTANCE, DDT.CYCLES_DISTANCE, DDT.FRAGMENTS_DISTANCE,
		  DDT.BREAKPOINT_DISTANCE, DDT.JACCARD_DISTANCE, DDT.TOTAL_COST]
colors = ["blue", "green", "red", "purple"]

def draw_total_cost(dict_metrics, outfile=None):
	# radial plot information
	r = []
	theta = []
	for key in dict_metrics[ht.DISTANCES][DDT.TOTAL_COST_DESCRIPTION]:
		theta.append(key)
		r.append(dict_metrics[ht.DISTANCES][key])

	df = pd.DataFrame(dict(
		r=r,
		theta=theta))

	fig = px.line_polar(df, r='r',
						theta='theta',
						line_close=True
						)
	# opacity = dict_metrics[ht.DISTANCES][DDT.TOTAL_COST] / (len(tokeep)-1)
	fig.update_traces(fill='toself', opacity=0.8)
	fig.update_layout(
		polar=dict(
			radialaxis=dict(
				visible=True,
				range=[0, 1]
			)),
		showlegend=True
	)

	if outfile:
		fig.write_image(outfile)
	else:
		# plot to stdin
		fig.show()


def draw_total_cost_table(dict_metrics, outfile=None):
	# table information
	table_val = []
	table_dist = []

	for key in dict_metrics[ht.DISTANCES]:
		if key in tokeep:
			table_dist.append(key)
			table_val.append(dict_metrics[ht.DISTANCES][key])

	fig = go.Figure(data=go.Table(
		header=dict(
			values=["Distance", "Value"],
			font=dict(color='darkslategray', size=14),
			line_color='gray',
			fill_color='white',
			align="left"
		),
		cells=dict(
			values=[table_dist, table_val],
			align="left",
			line_color='gray',
			fill_color='white',
			font=dict(color='darkslategray', size=12)
		)
	)
	)

	if outfile:
		fig.write_image(outfile)
	else:
		# plot to stdin
		fig.show()


def draw_breakpoints_creative(chr, start, end):
	angles = linspace(0 * pi + start, 1 * pi, 100)
	xs = abs(start - end) * cos(angles) + start
	ys = abs(start - end) * sin(angles) + start
	plt.plot(xs, ys, color='green')


def draw_arc(x, y, arc_start, arc_end, resolution=100, scale=False):
	"""
	arc_start and arc_end can be between 0 and 2, were
	0 means 0 degrees
	1 means 180 degrees
	2 means 360 degress
	"""
	angles = linspace(arc_start * pi, arc_end * pi, resolution)
	# xs = abs((x - y) / 2) * cos(angles) + x - (x - y) / 2
	# ys = abs((x - y) / 2) * sin(angles)
	xs = abs((x - y) / 2) * cos(angles) + x - (x - y) / 2

	if scale == True:
		ys = abs((x - y) / 2) * sin(angles)
	else:
		ys = sin(angles)
	return xs, ys


def draw_breakpoints(chr1, start, chr2, end, chr_offsets, color, alpha,
					 linesize=3, dotsize=10, w=5, ax=None, title="", flipped=False, scale=False):
	arc_start = 0
	arc_end = 1
	if flipped:
		arc_start = 1
		arc_end = 2

	# structure are on the same chromosome and no information about breakpoint match
	if chr1 == chr2:

		# get arc dots
		xs, ys = draw_arc(start, end, arc_start, arc_end, resolution=100, scale=scale)
		# set the chromosome offset
		xs = xs + chr_offsets[chr1]

		# draw arc
		ax.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
		ax.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
		ax.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
		ax.title.set_text(title)
		ax.set_yticks([])
		# else:
		# 	plt.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
		# 	plt.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
		# 	plt.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
		# 	plt.title(title)
		# plt.xlabel("")


def break_cn(df):
	"""
	Break copy number into the double number of bins
	"""
	cols = [ht.CHR, ht.START, ht.END, ht.CN, ht.TRACK]
	df_new = pd.DataFrame(columns=cols)
	j = 0
	chr_old = -1

	for i in range(0, df.shape[0]):
		chr, start, end, cn, type = df.loc[i, cols].tolist()
		len = end - start

		if chr_old == -1 or chr != chr_old:
			if chr_old != -1:
				df_new.loc[j, :] = [chr_old,
									df_new.loc[j - 1, ht.END],
									df_new.loc[j - 1, ht.END] + 100,
									df_new.loc[j - 1, ht.CN],
									df_new.loc[j - 1, ht.TRACK]]
				j = j + 1
			chr_old = chr

		df_new.loc[j, :] = [chr, start, start + math.floor(len / 2), cn, type]
		j = j + 1
		df_new.loc[j, :] = [chr, start + math.floor(len / 2) + 1, end, cn, type]
		j = j + 1

	# last element
	df_new.loc[j, :] = [chr_old,
						df_new.loc[j - 1, ht.END],
						df_new.loc[j - 1, ht.END] + 100,
						df_new.loc[j - 1, ht.CN],
						df_new.loc[j - 1, ht.TRACK]]

	return df_new


def draw_cn(cv_profile_t, cv_profile_r, chrlist, width=20, height=5, outfile=None):
	sns.set_style("whitegrid")
	sns.set_context("paper")

	# sns.set(rc={'figure.figsize': (width, height)})
	cv_profile_t[ht.TRACK] = ht.S1
	cv_profile_r[ht.TRACK] = ht.S2
	c_new = break_cn(cv_profile_t).append(break_cn(cv_profile_r), ignore_index=True)

	tracks = c_new[ht.TRACK].drop_duplicates().tolist()

	ncols = len(chrlist)
	fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True)
	ax = None

	if ncols == 1:
		ax = axs

	for i, ci in enumerate(chrlist):
		dict_x = {}
		dict_y = {}
		for j, t in enumerate(tracks):
			dict_x[t] = c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.END].tolist()
			dict_y[t] = c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.CN].tolist()

			if ncols == 1:
				ax = axs
			else:
				ax = axs[i]

			# flipp y-axis to negative for S2
			if t == ht.S2:
				dict_y[t] = [-y for y in dict_y[t]]

			ax.plot(dict_x[t], dict_y[t],
					drawstyle='steps',
					label=ci,
					linewidth=3,
					color=colors[j],
					alpha=0.5)
			ax.fill_between(dict_x[t], dict_y[t], color=colors[j], step="pre", alpha=0.1)

			ax.set_xlabel(ci, fontsize=12)
			ax.set_ylabel("")
			ax.set_xticklabels(np.array(ax.get_xticks()).astype(int),
							   rotation=90,
							   fontsize=12)
			ax.set_yticklabels(np.array(ax.get_yticks()).astype(int),
							   fontsize=12)

			if i == 0:
				ax.set_ylabel("cn", fontsize=12)

		# plot difference between the 2 profiles
		x_diff = dict_x[ht.S1]
		y_diff = [(np.log2(dict_y[ht.S1][k] + 1) - np.log2(abs(dict_y[ht.S2][k])+1)) for k in range(0, len(dict_y[ht.S1]))]
		ax.plot(x_diff, y_diff,
				drawstyle='steps',
				label=ht.S1 + "/" + ht.S2,
				linewidth=3,
				color="black",
				alpha=1)

	lines = [plt.Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors[:2] + ["black"]]
	ax.legend(lines, tracks[:2] + ["log2(" + ht.S1 + "/" + ht.S2 + ")"], fontsize=12)

	fig.show()
	if outfile:
		plt.savefig(outfile, bbox_inches='tight', dpi=300)


# plt.figure(figsize=(width, height))
# plt.plot(x, y + 2, drawstyle='steps', label='steps (=steps-pre)')
#
#
#
# ncols = c_new[ht.CHR].drop_duplicates().shape[0]
# d = {'color': ['green', 'blue'], "ls": ["-", "-"]}
#
#
# g = sns.FacetGrid(c_new,
# 				  col=ht.CHR,
# 				  height=3,
# 				  col_wrap=ncols,
# 				  hue="track",
# 				  hue_kws=d,
# 				  sharey="col",
# 				  sharex=None,
# 				  )
# g = g.map(sns.lineplot,
# 		   ht.END,
# 		   plot_col,
# 		   drawstyle='steps-pre',
# 		   alpha=0.5,
# 		   linewidth=3)
#
# # rename correctly the chr
# g.set_titles(col_template="")
# for i,val in enumerate(c_new[ht.CHR].drop_duplicates().tolist()):
# 	g.axes[i].set_xlabel(val)
# 	g.axes[i].set_xticklabels(np.array(g.axes[i].get_xticks()).astype(int), rotation = 90)
# g.add_legend()
# g.fig.suptitle(plot_col)


def plot_breakpoints_comparison(br_t, br_r, chr_offsets, breakpoint_match, chrlist, title="", width=20, height=5, outfile=None):
	"""

	Arguments:
		br_t (pd.DataFrame):
		br_r (pd.DataFrame):
		chr_offsets (dict):
		breakpoint_match (list): List of matching breakpoints ids  [(bt1,br1),(bt2,br2)]

	Returns:
	"""
	matched_br_t = [a[0] for a in breakpoint_match]
	matched_br_r = [a[1] for a in breakpoint_match]

	br_t[ht.TRACK] = ht.S1
	br_r[ht.TRACK] = ht.S2

	ncols = len(chrlist)
	fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True)

	ax = None
	if ncols != 1:
		ax_dict = {chr:axs[i] for i,chr in enumerate(chrlist)}

	for arr in [br_t, br_r]:
		for index, row in arr.iterrows():

			c1 = row[ht.CHR1]
			c2 = row[ht.CHR2]
			p1 = row[ht.START]
			p2 = row[ht.END]

			if c1 not in chrlist or c2 not in chrlist or c1 != c2:
				continue

			if ncols == 1:
				ax = axs
			else:
				ax = ax_dict[c1]

			# check if the breakpoint is a match
			color = PROPS.UNMATCHED[PROPS.COLOR]
			alpha = PROPS.UNMATCHED[PROPS.ALPHA]

			# set if flipped or not
			flipped = False
			if row[ht.TRACK] == ht.S2:
				flipped = True

			# if s1 and matched breakpoint
			if row[ht.TRACK] == ht.S1 and index in matched_br_t:
				color = colors[0]
				alpha = PROPS.MATCHED[PROPS.ALPHA]
				print(index)

			if row[ht.TRACK] == ht.S2 and index in matched_br_r:
				color = colors[1]
				alpha = PROPS.MATCHED[PROPS.ALPHA]

			draw_breakpoints(c1, p1, c2, p2, chr_offsets,
							color=color, alpha=alpha, ax=ax, flipped=flipped)



	ax.axhline(0, color="gray", alpha=0.5, ls="--")

