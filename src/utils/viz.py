import math
import pandas as pd
import numpy as np
from numpy import sin, cos, pi, linspace

import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

# custom module
from src.utils.utils import DDT
from src.utils.utils import HEADER as ht
from src.utils.utils import PROPS
from src.metrics import features


def draw_total_cost(dict_metrics, outfile=None):
	r = []
	theta = []
	for key in dict_metrics[ht.DISTANCES][DDT.TOTAL_COST_DESCRIPTION]:
		theta.append(key)
		r.append(dict_metrics[ht.DISTANCES][key])
	df = pd.DataFrame(dict(
		r=r,
		theta=theta))

	fig = px.line_polar(df, r='r', theta='theta', line_close=True)
	fig.update_traces(fill='toself')
	fig.update_layout(
		polar=dict(
			radialaxis=dict(
				visible=True,
				range=[0, 1]
			)),
		showlegend=False
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


def draw_breakpoints(chr1, start, chr2, end, chr_offsets, idx, breakpoints_list, color='green',
					 alpha=0.5, linesize=3, dotsize=10, w=5, ax=None, title="", flipped=False, scale=False):
	arc_start = 0
	arc_end = 1
	if flipped:
		arc_start = 1
		arc_end = 2

	# structure are on the same chromosome and no information about breakpoint match
	if chr1 == chr2:

		# breakpoints list set
		if breakpoints_list and idx not in breakpoints_list:
			# color unmatch gray
			color = PROPS.UNMATCHED[PROPS.COLOR]
			alpha = PROPS.UNMATCHED[PROPS.ALPHA]
		if not breakpoints_list or (breakpoints_list and len(breakpoints_list) == 0):
			# color unmatch gray
			color = PROPS.UNMATCHED[PROPS.COLOR]
			alpha = PROPS.UNMATCHED[PROPS.ALPHA]

		# get arc dots
		xs, ys = draw_arc(start, end, arc_start, arc_end, resolution=100, scale=scale)
		# set the chromosome offset
		xs = xs + chr_offsets[chr1]

		# draw arc
		if ax:
			ax.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
			ax.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
			ax.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
			ax.title.set_text(title)
			ax.set_yticks([])
		else:
			plt.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
			plt.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
			plt.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
			plt.title(title)
		plt.xlabel("")


def break_cn(df):
	"""
	Break copy number into the double number of bins
	"""
	cols = [ht.CHR, ht.START, ht.END, ht.CN, "track"]
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
									df_new.loc[j - 1, "track"]]
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
						df_new.loc[j - 1, "track"]]

	return df_new


def draw_cn(cv_profile_t, cv_profile_r, plot_col=ht.CN, width=20, height=5):
	sns.set_style("whitegrid")
	sns.set_context("paper")

	# sns.set(rc={'figure.figsize': (width, height)})
	cv_profile_t["track"] = "true"
	cv_profile_r["track"] = "reconstructed"
	c_new = break_cn(cv_profile_t).append(break_cn(cv_profile_r), ignore_index=True)

	chrlist = c_new[ht.CHR].drop_duplicates().tolist()
	colors = ["blue", "green", "red", "purple"]
	tracks = c_new["track"].drop_duplicates().tolist()

	ncols = len(chrlist)
	fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True)

	for i, ci in enumerate(chrlist):
		for j, t in enumerate(tracks):
			x = c_new.loc[(c_new[ht.CHR] == ci) & (c_new["track"] == t), ht.END].tolist()
			y = c_new.loc[(c_new[ht.CHR] == ci) & (c_new["track"] == t), ht.CN].tolist()
			axs[i].plot(x, y,
						drawstyle='steps',
						label=ci,
						linewidth=6,
						color=colors[j],
						alpha=0.7)
			axs[i].fill_between(x, y, color=colors[j], step="pre", alpha=0.2)

			axs[i].set_xlabel(ci, fontsize=12)
			axs[i].set_ylabel("")

			if i == 0:
				axs[i].set_ylabel("cn", fontsize=12)
			axs[i].set_xticklabels(np.array(axs[i].get_xticks()).astype(int),
								   rotation=90,
								   fontsize=12)
			axs[i].set_yticklabels(np.array(axs[i].get_yticks()).astype(int),
								   fontsize=12)
	axs[i].legend()


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


def plot_breakpoints_comparison(br_t, br_r, chr_offsets, breakpoint_match=None, ax=None, title="", scale=True):
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

	for index, row in br_t.iterrows():
		c1 = row["chr1"]
		c2 = row["chr2"]
		p1 = row["start"]
		p2 = row["end"]
		draw_breakpoints(c1, p1, c2, p2, chr_offsets, index, matched_br_t, color='green', alpha=0.5, ax=ax, scale=scale)

	for index, row in br_r.iterrows():
		c1 = row["chr1"]
		c2 = row["chr2"]
		p1 = row["start"]
		p2 = row["end"]
		draw_breakpoints(c1, p1, c2, p2, chr_offsets, index, matched_br_r, color='blue', alpha=0.5, ax=ax, title=title,
						 flipped=True, scale=scale)

	if ax:
		ax.axhline(0, color="gray", alpha=0.5, ls="--")
	else:
		plt.axhline(0, color="gray", alpha=0.5, ls="--")
