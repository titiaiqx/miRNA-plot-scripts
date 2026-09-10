# -*- coding: utf-8 -*-
"""
miRNA碱基组成
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import math

FONT_FAMILY = "Arial"
FONT_WEIGHT = "normal"

BASE_ORDER = ["A", "T", "C", "G"]
COLOR_A = "#4C9AC9"
COLOR_T = "#a6c9df"
COLOR_C = "#E2745E"
COLOR_G = "#F5B99E"
BASE_COLORS = {
    "A": COLOR_A,
    "T": COLOR_T,
    "C": COLOR_C,
    "G": COLOR_G
}

COMBINE_FIGSIZE = (12, 8)
SINGLE_FIGSIZE = (10, 6)
LEGEND_FIGSIZE = (6, 4)

BAR_WIDTH = 0.04
BASE_GAP = 0
GROUP_GAP = 0.04
BAR_BORDER_COLOR = "black"
BAR_BORDER_WIDTH = 1.5       

ANNOT_FONTSIZE = 20
ANNOT_FONTWEIGHT = "normal"
ANNOT_COLOR = "#000000"
ANNOT_OFFSET = 0.8

Y_AXIS_STEP = 5
Y_AXIS_MIN_START = 15
AXIS_LINE_WIDTH = 2
TICK_FONTSIZE = 24
TICK_FONTWEIGHT = "normal"

LABEL_FONTSIZE = 25
LABEL_FONTWEIGHT = "normal"

TITLE_FONTSIZE = 18
TITLE_FONTWEIGHT = "normal"
TITLE_PAD = 36
TITLE_LOC = "center"

LEGEND_FONTSIZE = 14
LEGEND_FONTWEIGHT = "normal"
LEGEND_BORDER_WIDTH = 1.5
LEGEND_FRAME_ON = True
LEGEND_MARKERSCALE = 1.2
LEGEND_BBOX_TO_anchor = (1.02, 1)

GRID_SHOW = True
GRID_COLOR = "#EAEAEA"
GRID_WIDTH = 1.2

SAVE_DPI = 300
# =============================================================================

def get_args():
    parser = argparse.ArgumentParser(description="miRNA ATCG总含量柱状图")
    parser.add_argument("-i", "--input", required=True, help="tsv文件")
    return parser.parse_args()

def load_miRNA_seq(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)

    def get_seqs(col):
        seqs = df.iloc[:, col].dropna().astype(str).str.strip().str.upper()
        seqs = [s for s in seqs if s]
        return seqs

    pre_seqs = get_seqs(11)
    mature_seqs = get_seqs(18)
    star_seqs = get_seqs(20)

    print(f"Pre-miRNA: {len(pre_seqs)}")
    print(f"Mature miRNA: {len(mature_seqs)}")
    print(f"Star miRNA: {len(star_seqs)}")
    return pre_seqs, mature_seqs, star_seqs

def calc_atcg_total_percent(seqs):
    A = T = C = G = 0
    for seq in seqs:
        A += seq.count("A")
        T += seq.count("T")
        C += seq.count("C")
        G += seq.count("G")
    total = A + T + C + G
    if total == 0:
        return [0.0, 0.0, 0.0, 0.0]
    pA = round(A / total * 100, 1)
    pT = round(T / total * 100, 1)
    pC = round(C / total * 100, 1)
    pG = round(G / total * 100, 1)
    return [pA, pT, pC, pG]

def auto_y_limit(all_pcts):
    max_val = max(all_pcts)
    y_max = math.ceil(max_val / Y_AXIS_STEP) * Y_AXIS_STEP
    y_min = Y_AXIS_MIN_START
    return y_min, y_max

def plot_combined(pre_pct, mat_pct, star_pct):
    plt.rcParams["font.family"] = FONT_FAMILY
    fig, ax = plt.subplots(figsize=COMBINE_FIGSIZE)
    n = len(BASE_ORDER)
    unit = BAR_WIDTH + BASE_GAP

    x_pre = np.arange(n) * unit
    x_mat = x_pre + n * unit + GROUP_GAP
    x_sta = x_mat + n * unit + GROUP_GAP

    for i, b in enumerate(BASE_ORDER):
        # Pre
        ax.bar(x_pre[i], pre_pct[i], BAR_WIDTH,
               color=BASE_COLORS[b], edgecolor=BAR_BORDER_COLOR, linewidth=BAR_BORDER_WIDTH)
        ax.text(x_pre[i], pre_pct[i] + ANNOT_OFFSET, f"{pre_pct[i]:.1f}",
                ha="center", va="bottom", fontsize=ANNOT_FONTSIZE,
                fontweight=ANNOT_FONTWEIGHT, color=ANNOT_COLOR)
        # Mature
        ax.bar(x_mat[i], mat_pct[i], BAR_WIDTH,
               color=BASE_COLORS[b], edgecolor=BAR_BORDER_COLOR, linewidth=BAR_BORDER_WIDTH)
        ax.text(x_mat[i], mat_pct[i] + ANNOT_OFFSET, f"{mat_pct[i]:.1f}",
                ha="center", va="bottom", fontsize=ANNOT_FONTSIZE,
                fontweight=ANNOT_FONTWEIGHT, color=ANNOT_COLOR)
        # Star
        ax.bar(x_sta[i], star_pct[i], BAR_WIDTH,
               color=BASE_COLORS[b], edgecolor=BAR_BORDER_COLOR, linewidth=BAR_BORDER_WIDTH)
        ax.text(x_sta[i], star_pct[i] + ANNOT_OFFSET, f"{star_pct[i]:.1f}",
                ha="center", va="bottom", fontsize=ANNOT_FONTSIZE,
                fontweight=ANNOT_FONTWEIGHT, color=ANNOT_COLOR)

    group_centers = [x_pre.mean(), x_mat.mean(), x_sta.mean()]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(["Pre-miRNA", "Mature miRNA", "Star miRNA"],
                       fontsize=TICK_FONTSIZE, fontweight=TICK_FONTWEIGHT)

    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE, width=AXIS_LINE_WIDTH)

    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in ["left", "bottom"])
        ax.spines[spine].set_linewidth(AXIS_LINE_WIDTH)

    if GRID_SHOW:
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=GRID_WIDTH)
        ax.set_axisbelow(True)

    ax.set_xlabel("miRNA Type", fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT)
    ax.set_ylabel("Base Percentage (%)", fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT)
    ax.set_title("ATCG Base Content Distribution of miRNA Sequences",
                 fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT,
                 pad=TITLE_PAD, loc=TITLE_LOC)

    all_vals = pre_pct + mat_pct + star_pct
    ax.set_ylim(auto_y_limit(all_vals))

    handles = [plt.Rectangle((0,0),1,1, fc=BASE_COLORS[b],
                             ec=BAR_BORDER_COLOR, lw=BAR_BORDER_WIDTH) for b in BASE_ORDER]
    leg = ax.legend(handles, BASE_ORDER,
              fontsize=LEGEND_FONTSIZE,
              frameon=LEGEND_FRAME_ON,
              markerscale=LEGEND_MARKERSCALE,
              loc="upper left",
              bbox_to_anchor=LEGEND_BBOX_TO_anchor)
    if leg.get_frame():
        leg.get_frame().set_linewidth(LEGEND_BORDER_WIDTH)

    plt.tight_layout()
    plt.savefig("miRNA_ATCG_content_combined.svg", format="svg", dpi=SAVE_DPI, bbox_inches="tight")
    plt.close()

def plot_single_figure(pct, seq_name):
    plt.rcParams["font.family"] = FONT_FAMILY
    fig, ax = plt.subplots(figsize=SINGLE_FIGSIZE)

    ax.bar(BASE_ORDER, pct,
           color=[BASE_COLORS[b] for b in BASE_ORDER],
           edgecolor=BAR_BORDER_COLOR,
           linewidth=BAR_BORDER_WIDTH)

    for i, val in enumerate(pct):
        ax.text(i, val + ANNOT_OFFSET, f"{val:.1f}",
                ha="center", va="bottom",
                fontsize=ANNOT_FONTSIZE,
                fontweight=ANNOT_FONTWEIGHT,
                color=ANNOT_COLOR)

    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE, width=AXIS_LINE_WIDTH)
    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in ["left", "bottom"])
        ax.spines[spine].set_linewidth(AXIS_LINE_WIDTH)

    if GRID_SHOW:
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=GRID_WIDTH)
        ax.set_axisbelow(True)

    ax.set_xlabel("Base Type", fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT)
    ax.set_ylabel("Base Percentage (%)", fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT)
    ax.set_title(f"{seq_name} ATCG Base Content",
                 fontsize=TITLE_FONTSIZE, fontweight=TITLE_FONTWEIGHT,
                 pad=TITLE_PAD, loc=TITLE_LOC)

    ax.set_ylim(auto_y_limit(pct))
    filename = f"miRNA_ATCG_content_{seq_name.lower()}.svg"
    plt.savefig(filename, format="svg", dpi=SAVE_DPI, bbox_inches="tight")
    plt.close()

def plot_alone_legend():
    fig, ax = plt.subplots(figsize=LEGEND_FIGSIZE)
    ax.axis("off")
    handles = [plt.Rectangle((0,0),1,1, fc=BASE_COLORS[b],
                             ec=BAR_BORDER_COLOR,
                             lw=BAR_BORDER_WIDTH) for b in BASE_ORDER]
    leg = ax.legend(handles, BASE_ORDER, fontsize=LEGEND_FONTSIZE,
                    frameon=LEGEND_FRAME_ON, loc="center")
    if leg.get_frame():
        leg.get_frame().set_linewidth(LEGEND_BORDER_WIDTH)
    plt.savefig("miRNA_ATCG_content_legend.svg", format="svg", dpi=SAVE_DPI, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    args = get_args()
    pre_seqs, mature_seqs, star_seqs = load_miRNA_seq(args.input)

    pre_pct = calc_atcg_total_percent(pre_seqs)
    mat_pct = calc_atcg_total_percent(mature_seqs)
    star_pct = calc_atcg_total_percent(star_seqs)

    plot_combined(pre_pct, mat_pct, star_pct)
    plot_single_figure(pre_pct, "Pre")
    plot_single_figure(mat_pct, "Mature")
    plot_single_figure(star_pct, "Star")
    plot_alone_legend()

    print("down!")