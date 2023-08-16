from numpy import sin, cos, pi, linspace
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

# custom module
from utils.utils import DDT
from utils.utils import HEADER as ht
from utils.utils import PROPS


def draw_total_cost(dict_metrics, outfile):
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

	# fig.show()
	fig.write_image(outfile)


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


def vizualize_cn(cv_profile_t, cv_profile_r, plot_col="coverage_mean", width=20, height=5):
	sns.set_style("whitegrid")
	sns.set_context("paper")
	# sns.set(rc={'figure.figsize': (width, height)})

	cv_profile_t["track"] = "true"
	cv_profile_r["track"] = "reconstructed"
	c = cv_profile_t.append(cv_profile_r)

	ncols = c["#chr"].drop_duplicates().shape[0]
	d = {'color': ['green', 'blue'], "ls": ["-", "-"]}

	g = sns.FacetGrid(c,
					  col=ht.CHR,
					  height=3,
					  col_wrap=ncols,
					  hue="track",
					  hue_kws=d,
					  sharey="col",
					  sharex=None
					  )
	ax = g.map(sns.lineplot,
			   ht.END,
			   plot_col,
			   drawstyle='steps-pre',
			   alpha=0.8,
			   linewidth=1.5)

	g.add_legend()
	g.fig.suptitle(plot_col)


def get_chromosome_offset(df_t, df_r):
	"""
	Get chromosome offsets for visualization purpose.

	Arguments:
		df_t (pd.DataFrame):
		df_r (pd.DataFrame:

	Returns:
	"""
	df = df_t.append(df_r, ignore_index=True)
	d1 = df[["#chr", "start"]].groupby(["#chr"]).max().reset_index()
	d1.columns = ["#chr", "pos"]

	d2 = df[["#chr", "end"]].groupby(["#chr"]).max().reset_index()
	d2.columns = ["#chr", "pos"]

	d3 = d1.append(d2, ignore_index=True).groupby(["#chr"]).max().reset_index()
	d3.columns = ["#chr", "pos"]

	chridx = d3["#chr"].drop_duplicates().tolist()

	offsets = {chridx[0]: 0}
	if len(chridx) == 1:
		return offsets

	l = d3.shape[0]
	cumm_offset = 0
	for i in range(1, l):
		max_pos = d3[d3["#chr"] == chridx[i - 1]].iloc[0, 1]
		cumm_offset += max_pos + 5000
		offsets[chridx[i]] = cumm_offset

	return offsets


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
