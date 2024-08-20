import math
import pprint

import pandas as pd
import numpy as np
from numpy import sin, cos, pi, linspace

import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.patches import ConnectionPatch, Circle

# custom module
from AmpliconComparison.utils.utils import DDT
from AmpliconComparison.utils.utils import HEADER as ht
from AmpliconComparison.utils.utils import PROPS
from AmpliconComparison.utils.utils import get_weight_distance, get_value_distance

tokeep = [DDT.HAMMING_NORM, DDT.COSINE_DISTANCE, DDT.CYCLES_DISTANCE, DDT.FRAGMENTS_DISTANCE,
		  DDT.BREAKPOINT_DISTANCE, DDT.JACCARD_DISTANCE, DDT.TOTAL_COST, DDT.COPYNUMBER_JC]
colors = ["#447604", "#3838CC", "red", "purple"]
custom_gray = "#D8DBE2"

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
		showlegend=True,
		font=dict(color='darkslategray', size=10)
	)

	if outfile:
		fig.write_image(outfile, scale=10, width=500, height=500)
		fig.write_image(outfile + ".svg", scale=10, width=500, height=500, format="svg")
	else:
		# plot to stdin
		fig.show()


def draw_total_cost_table(dict_metrics, outfile=None):
	# table information
	table_val = []
	table_dist = []
	table_val_w = []

	# pprint.pprint(dict_metrics)

	for d in dict_metrics[ht.DISTANCES]:
		if d in dict_metrics[ht.DISTANCES][DDT.TOTAL_COST_DESCRIPTION]:
			val = dict_metrics[ht.DISTANCES][d]
			weight = dict_metrics[ht.DISTANCES][DDT.TOTAL_COST_DESCRIPTION][d]
			table_dist.append(d)
			table_val.append(round(val, 2))
			table_val_w.append(round(val * weight, 2))

	# add total cost
	table_dist.append(DDT.TOTAL_COST)
	table_val.append("")
	table_val_w.append(dict_metrics[ht.DISTANCES][DDT.TOTAL_COST])

	fig = go.Figure(data=go.Table(
		header=dict(
			values=["Distance", "Value", "Weighted value"],
			font=dict(color='darkslategray', size=14),
			line_color=custom_gray,
			fill_color='white',
			align="left"
		),
		cells=dict(
			values=[table_dist, table_val, table_val_w],
			align="left",
			line_color=custom_gray,
			fill_color='white',
			font=dict(color='darkslategray', size=12)
		)
	)
	)

	if outfile:
		fig.write_image(outfile, scale=10, width=500, height=500)
		fig.write_image(outfile + ".svg", scale=10, width=500, height=500, format="svg")
	else:
		# plot to stdin
		fig.show()


def draw_breakpoints_creative(chr, start, end):
	angles = linspace(0 * pi + start, 1 * pi, 100)
	xs = abs(start - end) * cos(angles) + start
	ys = abs(start - end) * sin(angles) + start
	plt.plot(xs, ys, color='green')


def draw_arc(x, y, arc_start, arc_end, scale, max_value, resolution=100):
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

	if scale:
		if max_value:
			ymax = max_value - 0.1 * max_value
			# print("ymax",ymax)
			ys = min(max(20, abs((x - y) / 2)), ymax) * sin(angles)
		else:
			ys = min(20, abs((x - y) / 2)) * sin(angles)
	else:
		ys = sin(angles) + 20
	return xs, ys


def draw_breakpoints_cross_ax(chr1, start, chr2, end, max_value, scale, color, alpha,
							  ax1, ax2, fig, title="", linesize=3, dotsize=10, w=5, flipped=False):
	arc_start = 0
	arc_end = 1
	connstyle = 'arc3, rad=-0.5'

	if flipped:
		arc_start = 1
		arc_end = 2
		connstyle = 'arc3, rad=0.5'

	# xs, ys = draw_arc(start, end, arc_start, arc_end, scale, max_value, resolution=100)
	xy = (start, end)
	con1 = ConnectionPatch(xyA=(xy[0], 0),
						   coordsA=ax1.transData,
						   xyB=(xy[1], 0),
						   coordsB=ax2.transData,
						   color=color,
						   alpha=alpha,
						   linestyle="-",
						   linewidth=linesize,
						   connectionstyle=connstyle,
						   )
	fig.add_artist(con1)


# structure are on the same chromosome and no information about breakpoint match


# # get arc dots
# xs, ys = draw_arc(start, end, arc_start, arc_end, scale, max_value, resolution=100)
# # set the chromosome offset
# xs = xs + chr_offsets[chr1]
#
# if ax is not None:
# 	# draw arc
# 	ax.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
# 	ax.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
# 	ax.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
# 	ax.set_xlabel(chr1, fontsize=12)
# 	ax.title.set_text(title)
# 	# ax.set_yticks([])
# 	ax.set_xticklabels(np.array(ax.get_xticks()).astype(int),
# 					   rotation=90,
# 					   fontsize=12)
# 	if max_value:
# 		ax.set_ylim(-max_value, max_value)
# else:
# 	plt.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
# 	plt.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
# 	plt.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
# 	plt.title(title)
# 	plt.xlabel("")

def add_text(ax, fig, s1="", s2=""):
	trans = mtransforms.ScaledTranslation(10 / 72, -5 / 72, fig.dpi_scale_trans)
	ax.text(0.0, 1.0, ht.S1, transform=ax.transAxes + trans, fontsize='x-large', verticalalignment='top', zorder=10)
	trans = mtransforms.ScaledTranslation(10 / 72, 15 / 72, fig.dpi_scale_trans)
	ax.text(0.0, 0.0, ht.S2, transform=ax.transAxes + trans, fontsize='x-large', verticalalignment='top', zorder=10)

	ax.text(0.0, -0.3, ht.S1 + "=" + s1, transform=ax.transAxes + trans, fontsize='x-large', verticalalignment='top', zorder=10)
	ax.text(0.0, -0.35, ht.S2 + "=" + s2, transform=ax.transAxes + trans, fontsize='x-large', verticalalignment='top', zorder=10)

def draw_breakpoints(chr1, start, chr2, end, max_value, scale, color, alpha,
					 linesize=3, dotsize=10, w=5, ax=None, title="", flipped=False):
	arc_start = 0
	arc_end = 1
	if flipped:
		arc_start = 1
		arc_end = 2

	# structure are on the same chromosome and no information about breakpoint match
	if chr1 == chr2:

		# get arc dots
		xs, ys = draw_arc(start, end, arc_start, arc_end, scale, max_value, resolution=100)
		# set the chromosome offset
		# xs = xs + chr_offsets[chr1]

		if ax is not None:
			# draw arc
			ax.plot(xs, ys, color=color, alpha=alpha, lw=linesize)
			ax.scatter(xs[0], 0, s=w * dotsize, color=color, alpha=alpha)
			ax.scatter(xs[-1], 0, s=w * dotsize, color=color, alpha=alpha)
			ax.set_xlabel(chr1, fontsize=12)
			ax.title.set_text(title)
			# ax.set_yticks([])
			ax.set_xticklabels(np.array(ax.get_xticks()).astype(int),
							   rotation=90,
							   fontsize=12)
			if max_value:
				ax.set_ylim(-max_value, max_value)
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
									df_new.loc[j - 1, ht.END] + 1,
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
						df_new.loc[j - 1, ht.END] + 1,
						df_new.loc[j - 1, ht.CN],
						df_new.loc[j - 1, ht.TRACK]]

	return df_new


def draw_cn(cv_profile_t, cv_profile_r, chrlist, width=30, height=3, outfile=None, fig=None, axs=None, s1="", s2=""):
	sns.set_style("whitegrid")
	sns.set_context("paper")

	# sns.set(rc={'figure.figsize': (width, height)})
	cv_profile_t[ht.TRACK] = ht.S1
	cv_profile_r[ht.TRACK] = ht.S2
	c_new = pd.concat([break_cn(cv_profile_t),break_cn(cv_profile_r)], ignore_index=True)

	tracks = [ht.S1, ht.S2]
	ncols = len(chrlist)
	if fig is None or axs is None:
		points_new = concat_all_breakpoints(cv_profile_t, cv_profile_r, ht.CHR, ht.CHR, ht.CHR, ht.START, ht.END)
		compute_ratios = plot_ratios(chrlist, points_new)
		fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True, gridspec_kw={'width_ratios': compute_ratios})
	ax = None

	if ncols == 1:
		ax = axs

	for i, ci in enumerate(chrlist):
		dict_x = {}
		dict_y = {}
		for j, t in enumerate(tracks):
			dict_x[t] = c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.END].tolist()
			dict_y[t] = c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.CN].tolist()

			dict_x[t] = [c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.START].tolist()[0]] + dict_x[t]
			dict_y[t] = [c_new.loc[(c_new[ht.CHR] == ci) & (c_new[ht.TRACK] == t), ht.CN].tolist()[0]] + dict_y[t]

			if ncols == 1:
				ax = axs
			else:
				ax = axs[i]

			# flipp y-axis to negative for S2
			color = colors[0]
			if t == ht.S2:
				dict_y[t] = [-y for y in dict_y[t]]
				color = colors[1]

			ax.plot(dict_x[t], dict_y[t],
					drawstyle='steps',
					label=ci,
					linewidth=3,
					color="gray",
					alpha=0.8)
			ax.fill_between(dict_x[t], dict_y[t], color=custom_gray, step="pre", alpha=0.8)

			ax.set_xlabel(ci, fontsize=12)
			ax.set_ylabel("")
			ax.set_xticklabels(np.array(ax.get_xticks()).astype(int),
							   rotation=90,
							   fontsize=12)
			ax.set_yticklabels(np.array(ax.get_yticks()).astype(int),
							   fontsize=12)

			if i == 0:
				ax.set_ylabel("estimated copy-number", fontsize=12)
				add_text(ax, fig, s1=s1, s2=s2)

		# plot difference between the 2 profiles
		x_diff = dict_x[ht.S1]
		y_diff = [(np.log2(dict_y[ht.S1][k] + 0.01) - np.log2(abs(dict_y[ht.S2][k]) + 0.01)) for k in
				  range(0, len(dict_y[ht.S1]))]

		ax.axhline(0, color="red", alpha=1, ls="--")
		ax.plot(x_diff, y_diff,
				drawstyle='steps',
				label=ht.S1 + "/" + ht.S2,
				linewidth=4,
				color="black",
				alpha=0.8,
		        zorder=1)
	if outfile:
		fig.savefig(outfile, bbox_inches='tight', dpi=600)
		fig.savefig(outfile+".svg", bbox_inches='tight', format="svg")
	else:
		fig.show()


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

def concat_all_breakpoints(df1, df2, achr1, achr2, achr, astart, aend):
	br_t_1 = pd.DataFrame(np.concatenate((df1[[achr1, astart]].values,
										  df1[[achr2, aend]].values),
										 axis=0))
	br_r_1 = pd.DataFrame(np.concatenate((df2[[achr1, astart]].values,
										  df2[[achr2, aend]].values),
										 axis=0))

	points_new = pd.concat([br_t_1,br_r_1], ignore_index=True).drop_duplicates()
	points_new.columns = [achr, astart]
	return points_new


def plot_ratios(chr_list, points_new, min_subplot_ratio=0.03):

	ratios = []

	for c in chr_list:
		df_gr = points_new[points_new[ht.CHR]==c].groupby(by=ht.CHR, observed=True,as_index=False).agg({ht.START: [np.min, np.max]})
		min_, max_ = df_gr.iloc[0,1], df_gr.iloc[0,2]
		ratios.append(max_ - min_)

	new_ratios = np.array(ratios) / np.array(ratios).sum()

	# make sure some very small bins do not go unter 0.001 ratio
	max_ratio = np.argmax(ratios)
	count_shrink = 0
	new_ratios_adj = []
	for n in new_ratios:
		if n > min_subplot_ratio:
			new_ratios_adj.append(n)
		else:
			new_ratios_adj.append(min_subplot_ratio)
			count_shrink += 1
	new_ratios_adj[max_ratio] = new_ratios_adj[max_ratio] - count_shrink * min_subplot_ratio

	return new_ratios_adj

def plot_breakpoints_location(br_t, br_r, max_y, chrlist, width=30, height=5, s1="", s2=""):
	"""
	Plot as dots all the breakpoints
	"""

	sns.set_style("whitegrid")
	sns.set_context("paper")

	points_new = concat_all_breakpoints(br_t, br_r, ht.CHR1, ht.CHR2, ht.CHR, ht.START, ht.END)
	ncols = len(chrlist)
	compute_ratios = plot_ratios(chrlist, points_new)
	fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True, gridspec_kw={'width_ratios': compute_ratios})

	for i, ci in enumerate(chrlist):

		x = points_new.loc[points_new[ht.CHR] == ci, ht.START].tolist()
		y = np.zeros(len(x))

		if ncols == 1:
			ax = axs
		else:
			ax = axs[i]

		ax.scatter(x, y,
				   label=ci,
				   color=custom_gray,
				   alpha=0.5)

		ax.set_xlabel(ci, fontsize=12)
		ax.set_ylim((-max_y, max_y))
		ax.set_ylabel("")
		ax.set_xticklabels(np.array(ax.get_xticks()).astype(int),
						   rotation=90,
						   fontsize=12)
		ax.set_yticklabels(np.array(ax.get_yticks()).astype(int),
						   fontsize=12)

		# label samples
		if i == 0:
			add_text(ax, fig, s1=s1, s2=s2)
	return fig, axs


def plot_breakpoints_comparison(br_t, br_r, breakpoint_match, chrlist,
								width=30, height=3, max_value=None, scale=True, fig=None, axs=None, s1="",s2=""):
	"""

	Arguments:
		br_t (pd.DataFrame):
		br_r (pd.DataFrame):
		breakpoint_match (list): List of matching breakpoints ids  [(bt1,br1),(bt2,br2)]
		chrlist:
	Returns:
	"""
	matched_br_t = [a[0] for a in breakpoint_match]
	matched_br_r = [a[1] for a in breakpoint_match]

	br_t[ht.TRACK] = ht.S1
	br_r[ht.TRACK] = ht.S2

	ax = None
	ncols = len(chrlist)
	if fig is None and axs is None:
		points_new = concat_all_breakpoints(br_t, br_r, ht.CHR1, ht.CHR2, ht.CHR, ht.START, ht.END)
		compute_ratios = plot_ratios(chrlist, points_new)
		fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True, gridspec_kw={'width_ratios': compute_ratios})
		ax = axs[0] if ncols > 1 else axs
		ax.set_yticklabels([])
		add_text(ax, fig, s1=s2, s2=s2)

	if ncols > 1:
		ax_dict = {chr: axs[i] for i, chr in enumerate(chrlist)}

	for arr in [br_t, br_r]:
		for index, row in arr.iterrows():

			c1 = row[ht.CHR1]
			c2 = row[ht.CHR2]
			p1 = row[ht.START]
			p2 = row[ht.END]

			if c1 != c2 or c1 not in chrlist or c2 not in chrlist:
				continue

			# make sure the chr are ordered
			if chrlist.index(c1) > chrlist.index(c2):
				# swap
				c1 = row[ht.CHR2]
				c2 = row[ht.CHR1]
				p1 = row[ht.END]
				p2 = row[ht.START]

			if ncols == 1:
				ax = axs
			else:
				ax = ax_dict[c1]

			# check if the breakpoint is a match
			color = PROPS.UNMATCHED[PROPS.COLOR]
			alpha = PROPS.UNMATCHED[PROPS.ALPHA]

			flipped = True if row[ht.TRACK] == ht.S2 else False

			# if s1 and matched breakpoint
			if row[ht.TRACK] == ht.S1 and index in matched_br_t:
				color = colors[0]
				alpha = PROPS.MATCHED[PROPS.ALPHA]

			if row[ht.TRACK] == ht.S2 and index in matched_br_r:
				color = colors[1]
				alpha = PROPS.MATCHED[PROPS.ALPHA]

			draw_breakpoints(c1, p1, c2, p2, max_value, scale,
							 color, alpha, ax=ax, flipped=flipped)

			if max_value:
				ax.set_ylim(-max_value, max_value)

	# add links connecting chromosomes
	for arr in [br_t, br_r]:

		for index, row in arr.iterrows():

			c1 = row[ht.CHR1]
			c2 = row[ht.CHR2]
			p1 = row[ht.START]
			p2 = row[ht.END]

			# skip breakpoints from the same chromosome
			if c1 == c2:
				continue

			# make sure the chr are ordered
			ax1 = ax_dict[c1]
			ax2 = ax_dict[c2]
			if chrlist.index(c1) > chrlist.index(c2):
				# swap
				c1 = row[ht.CHR2]
				c2 = row[ht.CHR1]
				p1 = row[ht.END]
				p2 = row[ht.START]
				ax1 = ax_dict[c1]
				ax2 = ax_dict[c2]

			# check if the breakpoint is a match
			color = PROPS.UNMATCHED[PROPS.COLOR]
			alpha = PROPS.UNMATCHED[PROPS.ALPHA]

			# set if flipped or not
			flipped = True if row[ht.TRACK] == ht.S2 else False

			# if s1 and matched breakpoint
			if row[ht.TRACK] == ht.S1 and index in matched_br_t:
				color = colors[0]
				alpha = PROPS.MATCHED[PROPS.ALPHA]

			if row[ht.TRACK] == ht.S2 and index in matched_br_r:
				color = colors[1]
				alpha = PROPS.MATCHED[PROPS.ALPHA]

			draw_breakpoints_cross_ax(c1, p1, c2, p2, max_value, scale, color, alpha,
									  ax1, ax2, fig, flipped=flipped)


def plot_combined(br_t, br_r, cn_profile_t, cn_profile_r, breakpoint_matches, chrlist, max_coverage, outfile=None, s1="", s2=""):

	fig, axs = plot_breakpoints_location(br_t, br_r, max_coverage, chrlist)
	draw_cn(cn_profile_t, cn_profile_r, chrlist, fig=fig, axs=axs, s1=s1, s2=s2)
	plot_breakpoints_comparison(br_t, br_r, breakpoint_matches, chrlist, fig=fig, axs=axs, max_value=max_coverage, scale=True, s1=s1, s2=s2)
	if outfile:
		fig.savefig(outfile, bbox_inches='tight', dpi=600)
		fig.savefig(outfile + ".svg", bbox_inches='tight', format="svg")
	else:
		fig.show()
