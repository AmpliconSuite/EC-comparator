"""
threadsviz.py

Draws "reconstruction threads": the path/cycle structure of ecDNA (or other
rearranged) reconstructions, showing how genomic segments are joined together
in order, with strand-aware directionality, and support for both circular
(default) and linear (ISCYCLIC == False) structures.

Input format (tab-separated), one row per segment:
    Chromosome  Start    End circ_id  CN Stranded  iscyclic

Connections between consecutive segments are drawn as flat-topped, round-
cornered "bracket" links (a single cubic Bezier per link) rather than a
circular bulge, so link height no longer grows with genomic distance -
every link in a row sits at the same plateau height, and only curves where
it meets a segment.
"""

import io
import pprint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path

from eccomparator.utils.utils import HEADER

# ---- styling (reusing the same look/conventions as eccomparator's viz.py) ----
STRUCTURE_COLORS = ["#447604", "#3838CC", "#B23A48", "#8E44AD"]
CUSTOM_GRAY = "#D8DBE2"
ROW_HEIGHT = 0.35
BAR_HEIGHT = 0.22
PLATEAU_OFFSET = 0.22   # how far above/below a row's baseline the flat link plateau sits


def load_structure(path_or_buffer):
    """Read a tab-separated segment table into a DataFrame."""
    df = pd.read_csv(path_or_buffer, sep="\t")
    return df


def _parse_bool(val):
    """Parse a boolean that may arrive as an actual bool, or as the
    string 'True'/'False' (as read from a TSV column, e.g. via pandas)."""
    if isinstance(val, str):
        return val.strip().lower() in ("true", "1", "yes")
    return bool(val)


def _segment_paths(df):
    """
    Group segments by circ_id, preserving row order as the path order.
    Returns a list of dicts: {circ_id, iscyclic, segments: [row, ...]}
    """
    paths = []
    for circ_id, group in df.groupby("circ_id", sort=False):
        paths.append({
            HEADER.CIRC_ID: circ_id,
            HEADER.ISCYCLIC: _parse_bool(group[HEADER.ISCYCLIC].iloc[0]),
            HEADER.SEGMENTS: group.to_dict("records"),
        })
    return paths


def _chrom_order_and_ranges(dfs):
    """Every chromosome referenced across all given DataFrames, with its min/max span."""
    ranges = {}
    for df in dfs:
        for _, row in df.iterrows():
            c = row[HEADER.CHR]
            lo, hi = ranges.get(c, (row[HEADER.START], row[HEADER.END]))
            ranges[c] = (min(lo, row[HEADER.START]), max(hi, row[HEADER.END]))
    chrom_order = sorted(ranges, key=lambda c: (len(c), c))
    return chrom_order, ranges


def _width_ratios(chrom_order, ranges, min_ratio=0.10):
    spans = np.array([ranges[c][1] - ranges[c][0] for c in chrom_order], dtype=float)
    ratios = spans / spans.sum()
    ratios = np.maximum(ratios, min_ratio)
    ratios = ratios / ratios.sum()
    return ratios


def _entry_exit(seg):
    """
    Entry/exit coordinates of a segment, respecting strand.
    '+' strand: enters at start, exits at end   (traversed left -> right)
    '-' strand: enters at end,   exits at start  (traversed right -> left)
    """
    if seg[HEADER.STRAND] == "-":
        return seg[HEADER.END], seg[HEADER.START]
    return seg[HEADER.START], seg[HEADER.END]


def _make_axes_row(fig, gridspec_cells, chrom_order):
    axs = [fig.add_subplot(cell) for cell in gridspec_cells]
    return {c: axs[i] for i, c in enumerate(chrom_order)}


def _flat_bracket(fig, ax_a, x_a, y_a, ax_b, x_b, y_b, plateau_y, color,
                   linestyle="-", lw=1.8, alpha=0.8, zorder=2):
    """
    Draw a flat-topped, round-cornered link between (x_a, y_a) on ax_a and
    (x_b, y_b) on ax_b. Implemented as one cubic Bezier per pair, with control
    points pulled vertically to `plateau_y` - this gives a shape that is flat
    across the middle and rounds off smoothly right at the segment ends,
    regardless of how far apart x_a/x_b are.

    Works identically whether ax_a is ax_b (same chromosome) or not (cross-
    chromosome), since everything is expressed in figure-fraction coordinates.
    """
    to_fig = fig.transFigure.inverted()

    plateau_a = min(plateau_y + np.random.uniform(0.3, 0.8), y_a + 0.8)
    plateau_b = min(plateau_y + np.random.uniform(0.3, 0.8), y_b + 0.8)

    pA = to_fig.transform(ax_a.transData.transform((x_a, y_a)))
    pA_plateau = to_fig.transform(ax_a.transData.transform((x_a, plateau_a)))
    pB = to_fig.transform(ax_b.transData.transform((x_b, y_b)))
    pB_plateau = to_fig.transform(ax_b.transData.transform((x_b, plateau_b)))

    verts = [tuple(pA), tuple(pA_plateau), tuple(pB_plateau), tuple(pB)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = PathPatch(path, transform=fig.transFigure, facecolor="none",
                       edgecolor=color, linestyle=linestyle, lw=lw,
                       alpha=alpha, capstyle="round", zorder=zorder)
    fig.add_artist(patch)


EDGE_COLOR = "#666666"
DUP_LANE_HEIGHT = 0.25   # vertical separation between repeated visits to the same region


def _edge_linewidth(cn):
    """Edge thickness scales with copy number (same convention as pyGenomeTracks: 0.5*sqrt(score))."""
    return 0.5 + 0.5 * np.sqrt(max(cn, 0))


def _assign_lanes(segs):
    """
    Greedy first-fit lane assignment, processed in path order (not sorted by
    position): the first segment gets lane 0; each subsequent segment is
    placed in the lowest-numbered existing lane it does NOT overlap with
    (same chromosome and overlapping coordinates), or a new lane if none of
    the existing lanes are free. Two non-overlapping segments can share a
    lane, so the number of lanes used is kept minimal for the given order.

    Returns a list of lane indices, one per segment in `segs`.
    """
    lanes = []  # each entry: list of (x0, x1, chrom) already placed in that lane
    lane_of = []
    for seg in segs:
        x0, x1 = min(seg[HEADER.START], seg[HEADER.END]), max(seg[HEADER.START], seg[HEADER.END])
        chrom = seg[HEADER.CHR]
        placed_lane = None
        for lane_idx, occupied in enumerate(lanes):
            overlaps_lane = any(
                chrom == o_chrom and x0 < o_x1 and o_x0 < x1
                for o_x0, o_x1, o_chrom in occupied
            )
            if not overlaps_lane:
                placed_lane = lane_idx
                break
        if placed_lane is None:
            placed_lane = len(lanes)
            lanes.append([])
        lanes[placed_lane].append((x0, x1, chrom))
        lane_of.append(placed_lane)
    return lane_of


def _draw_structure_row(fig, ax_dict, df, color, sign=1, label_ax=None, name=""):
    """
    Draw one structure's segments + strand arrows + reconstruction links onto
    an existing dict of {chrom: ax}. `sign` controls whether rows stack
    upward (+1) or downward (-1) from that row's own y=0 baseline.

    Segments that overlap each other in genomic coordinates (including exact
    duplicates) are NOT drawn on top of one another: each is assigned to its
    own "lane" via a greedy first-fit algorithm, so any two segments sharing
    a lane are guaranteed not to overlap, and the number of lanes used is
    kept as small as possible for the given path order.
    """
    paths = _segment_paths(df)

    y_cursor = 0.0
    for row_idx, path in enumerate(paths):
        segs = path[HEADER.SEGMENTS]

        # assign each segment to the lowest-numbered lane it doesn't overlap
        # with; non-overlapping segments can share a lane, overlapping ones
        # (including exact duplicates) never do
        lane_of = _assign_lanes(segs)
        n_lanes = max(lane_of) + 1 if lane_of else 1

        y_cursor += ROW_HEIGHT
        base_y = sign * y_cursor
        seg_y = [base_y + sign * lane * DUP_LANE_HEIGHT for lane in lane_of]
        y_cursor += (n_lanes - 1) * DUP_LANE_HEIGHT  # reserve this row's extra lanes before the next row starts

        for seg, y in zip(segs, seg_y):
            ax = ax_dict[seg[HEADER.CHR]]
            x0, x1 = min(seg[HEADER.START], seg[HEADER.END]), max(seg[HEADER.START], seg[HEADER.END])
            rect = Rectangle((x0, y - BAR_HEIGHT / 2), x1 - x0, BAR_HEIGHT,
                              facecolor=color, edgecolor="none", alpha=0.9, zorder=3)
            ax.add_patch(rect)

            mid = (x0 + x1) / 2
            dx = (x1 - x0) * 0.12 if seg[HEADER.STRAND] == "+" else -(x1 - x0) * 0.12
            ax.annotate("", xy=(mid + dx, y), xytext=(mid - dx, y),
                        arrowprops=dict(arrowstyle="-|>", color="white",
                                         lw=1.4, shrinkA=0, shrinkB=0),
                        zorder=4)
        
        n = len(segs)
        pairs = [(segs[k], seg_y[k], segs[k + 1], seg_y[k + 1]) for k in range(n - 1)]
        if path[HEADER.ISCYCLIC]==True and n >= 1:
            pairs.append((segs[-1], seg_y[-1], segs[0], seg_y[0]))

        for a, y_a, b, y_b in pairs:
            _, exit_a = _entry_exit(a)
            entry_b, _ = _entry_exit(b)
            cross_chrom = a[HEADER.CHR] != b[HEADER.CHR]
            plateau_y = (y_a + y_b) / 2 + sign * PLATEAU_OFFSET
            cn_avg = (a[HEADER.CN] + b[HEADER.CN]) / 2
            _flat_bracket(fig, ax_dict[a[HEADER.CHR]], exit_a, y_a,
                           ax_dict[b[HEADER.CHR]], entry_b, y_b, plateau_y,
                           color=EDGE_COLOR, lw=_edge_linewidth(cn_avg),
                           linestyle="--" if cross_chrom else "-")

        if label_ax is not None:
            label = f"{name} \u00b7 circ {path[HEADER.CIRC_ID]}" + ("" if path[HEADER.ISCYCLIC] else " (linear)")
            label_ax.annotate(label, xy=(0, base_y), xycoords=label_ax.get_yaxis_transform(),
                               xytext=(-8, 0), textcoords="offset points",
                               ha="right", va="center", fontsize=9, clip_on=False)


def _style_chrom_axes(axs, chrom_order, ranges, y_range):
    lo_y, hi_y = y_range
    for i, c in enumerate(chrom_order):
        ax = axs[i]
        lo, hi = ranges[c]
        pad = 0.08 * (hi - lo if hi > lo else 1)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo_y, hi_y)
        ax.set_xlabel(c, fontsize=12)
        ax.axhline(0, color=CUSTOM_GRAY, lw=1, zorder=0)
        ax.set_yticks([])
        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)


def draw_reconstruction_threads(structures, outfile=None, width=16, height=6):
    """
    Standalone combined view: all given structures overlaid on one row of
    per-chromosome axes (first structure stacked upward, second downward,
    alternating further out for additional structures).

    structures: dict of {structure_name: DataFrame}, e.g. {"s1": df1, "s2": df2}
    """
    chrom_order, ranges = _chrom_order_and_ranges(structures.values())
    ratios = _width_ratios(chrom_order, ranges)
    ncols = len(chrom_order)

    fig, axs = plt.subplots(1, ncols, figsize=(width, height), sharey=True,
                             gridspec_kw={"width_ratios": ratios})
    if ncols == 1:
        axs = [axs]
    ax_dict = {c: axs[i] for i, c in enumerate(chrom_order)}

    band_sign = {name: (1 if idx % 2 == 0 else -1) for idx, name in enumerate(structures)}
    max_rows = max(len(_segment_paths(df)) for df in structures.values())
    ylim = (max_rows + 1) * ROW_HEIGHT * 2.2
    _style_chrom_axes(axs, chrom_order, ranges, (-ylim, ylim))

    legend_handles = []
    for s_idx, (name, df) in enumerate(structures.items()):
        color = STRUCTURE_COLORS[s_idx % len(STRUCTURE_COLORS)]
        _draw_structure_row(fig, ax_dict, df, color, sign=band_sign[name],
                             label_ax=axs[0], name=name)
        legend_handles.append(plt.Line2D([0], [0], color=color, lw=4, label=name))

    fig.legend(handles=legend_handles, loc="upper center", ncol=len(structures),
               frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Reconstruction threads", y=1.08, fontsize=14)
    fig.tight_layout()

    if outfile:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, axs


def draw_reconstruction_overview(structures, outfile=None, width=16, in_per_track=1.3,
                                  base_in=1.2, show_combined=False):
    """
    Multi-row figure, all columns (chromosomes) aligned:
      Row 0 (optional): the combined overlay (all structures together, as in
                         draw_reconstruction_threads) - only shown if
                         show_combined=True.
      One row per structure, alone, centered on y=0.

    Each row's physical height scales with how many  tracks it needs
    to stack, so bars/links/labels keep the same visual scale regardless of
    how many rows the overall figure has.

    structures: dict of {structure_name: DataFrame}, e.g. {"s1": df1, "s2": df2}
    """
    names = list(structures.keys())
    chrom_order, ranges = _chrom_order_and_ranges(structures.values())
    ratios = _width_ratios(chrom_order, ranges)
    ncols = len(chrom_order)

    n_rows_per_name = {name: len(_segment_paths(structures[name])) for name in names}
    combined_tracks = 2 * max(n_rows_per_name.values())          # stacks both up and down
    row_track_counts = ([combined_tracks] if show_combined else []) + \
                       [n_rows_per_name[name] for name in names]
    row_heights_in = [base_in + in_per_track * t for t in row_track_counts]
    total_height = sum(row_heights_in) + 1.4  # + headroom for suptitle/legend

    fig = plt.figure(figsize=(width, total_height))
    gs = fig.add_gridspec(len(row_heights_in), ncols, width_ratios=ratios,
                           height_ratios=row_heights_in, hspace=0.9,
                           top=1.0 - 1.4 / total_height)

    legend_handles = [plt.Line2D([0], [0], color=STRUCTURE_COLORS[i % len(STRUCTURE_COLORS)],
                                  lw=4, label=name) for i, name in enumerate(names)]

    row_offset = 0
    if show_combined:
        # ---- row 0: combined overlay ----
        row0 = _make_axes_row(fig, [gs[0, i] for i in range(ncols)], chrom_order)
        band_sign = {name: (1 if idx % 2 == 0 else -1) for idx, name in enumerate(names)}
        max_rows = max(n_rows_per_name.values())
        axs0 = [row0[c] for c in chrom_order]
        _style_chrom_axes(axs0, chrom_order, ranges, (-(max_rows + 1) * ROW_HEIGHT * 2.2,
                                                        (max_rows + 1) * ROW_HEIGHT * 2.2))
        for s_idx, name in enumerate(names):
            color = STRUCTURE_COLORS[s_idx % len(STRUCTURE_COLORS)]
            _draw_structure_row(fig, row0, structures[name], color,
                                 sign=band_sign[name], label_ax=row0[chrom_order[0]], name=name)
        axs0[0].set_title("Combined", loc="left", fontsize=11, style="italic", color="dimgray")
        row_offset = 1

    fig.legend(handles=legend_handles, loc="upper center", ncol=len(names),
               frameon=False, bbox_to_anchor=(0.5, 1.0 - 0.55 / total_height))

    # ---- one row per structure, alone ----
    for s_idx, name in enumerate(names):
        color = STRUCTURE_COLORS[s_idx % len(STRUCTURE_COLORS)]
        row = _make_axes_row(fig, [gs[row_offset + s_idx, i] for i in range(ncols)], chrom_order)
        n_rows = n_rows_per_name[name]
        axs_row = [row[c] for c in chrom_order]
        top = (n_rows + 0.6) * ROW_HEIGHT * 2.4
        bottom = -0.6 * ROW_HEIGHT * 2.4
        _style_chrom_axes(axs_row, chrom_order, ranges, (bottom, top))
        _draw_structure_row(fig, row, structures[name], color, sign=1,
                             label_ax=row[chrom_order[0]], name=name)
        axs_row[0].set_title(f"Reconstruction: {name}", loc="left", fontsize=11,
                              style="italic", color="dimgray")

    fig.suptitle("Reconstruction threads", y=1.0 - 0.05 / total_height, fontsize=15)

    if outfile:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig


if __name__ == "__main__":
    s1_tsv = """Chromosome\tStart\tEnd\tStrand\t\tCN\tiscyclic
chr1\t1000\t2000\t+\t1\t10\tTrue
chr1\t7000\t8000\t+\t1\t10\tTrue
chr1\t600\t2500\t-\t1\t10\tTrue
chr1\t1500\t2100\t+\t2\t20\tTrue
chr1\t9000\t11100\t-\t2\t20\tTrue
chr2\t100\t2000\t+\t2\t20\tTrue
"""
    s2_tsv = """Chromosome\tStart\tEnd\tStrand\t\tCN\tiscyclic
chr1\t500\t2500\t-\t1\t10\tTrue
chr1\t7200\t8050\t+\t1\t10\tTrue
chr1\t1500\t2100\t+\t2\t20\tTrue
chr1\t9000\t11100\t+\t2\t20\tTrue
chr2\t100\t2000\t+\t2\t20\tTrue
"""
    df_s1 = load_structure(io.StringIO(s1_tsv))
    df_s2 = load_structure(io.StringIO(s2_tsv))

    draw_reconstruction_overview({"s1": df_s1, "s2": df_s2},
                                  outfile="reconstruction_overview_demo.png")
