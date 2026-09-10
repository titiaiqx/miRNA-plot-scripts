# -*- coding: utf-8 -*-
"""
miRNA首碱基
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math

CUSTOM_STYLE = {
    "figure_dpi": 300,
    "font_family": "Arial",
    "base_order": ["A", "T", "C", "G"],
    "base_colors": {
        "A": "#4C9AC9",
        "T": "#a6c9df",
        "C": "#E2745E",
        "G": "#F5B99E"
    },

    # --------------- 画布尺寸 ---------------
    "combine_figsize": (12, 8),
    "single_figsize": (10, 6),
    "legend_figsize": (6, 4),

    "bar_width": 0.04,
    "base_gap": 0,
    "group_gap": 0.04,

    "bar_border_color": "black",
    "bar_border_width": 1.5,
    "bar_alpha": 1.0,

    "annot_fontsize": 20,
    "annot_fontweight": "normal",
    "annot_color": "#000000",
    "annot_offset": 0.8,

    "y_axis_extend_multiple": 5,
    "axis_line_width": 2,
    "axis_spines_show": ["left", "bottom"],

    "grid_show": True, 
    "grid_color": "#EAEAEA",
    "grid_width": 1.2,

    "tick_fontsize": 24,
    "tick_fontweight": "normal",

    "label_fontsize": 25,
    "label_fontweight": "normal",
    "y_label_fontsize": 28,

    "title_fontsize": 18,
    "title_fontweight": "normal",
    "title_pad": 36,
    "title_loc": "center",

    "legend_fontsize": 14,
    "legend_fontweight": "normal",
    "legend_border_width": 1.5,
    "legend_frame_on": True,
    "legend_markerscale": 1.2,
    "legend_loc": "upper left",
    "legend_bbox_to_anchor": (1.02, 1),

    "pre_col": 11,
    "mature_col": 18,
    "star_col": 20,
}

plt.rcParams['font.family'] = CUSTOM_STYLE["font_family"]
plt.rcParams['axes.unicode_minus'] = False

def parse_args():
    parser = argparse.ArgumentParser(description="miRNA首碱基")
    parser.add_argument("-i", "--input", required=True, help="tsv文件")
    return parser.parse_args()

def load_first_base(file_path, style):
    df = pd.read_csv(file_path, sep="\t", header=None)
    def get_first_bases(col_idx):
        seqs = df.iloc[:, col_idx].dropna().astype(str).str.strip()
        seqs = [s for s in seqs if s and len(s) > 0]
        first_bases = [s[0].upper() for s in seqs]
        return first_bases
    pre_bases = get_first_bases(style["pre_col"])
    mature_bases = get_first_bases(style["mature_col"])
    star_bases = get_first_bases(style["star_col"])
    return pre_bases, mature_bases, star_bases

def calc_base_percent(bases, base_order):
    if not bases:
        return [0.0]*len(base_order)
    total = len(bases)
    cnt = Counter(bases)
    return [round(cnt.get(b, 0)/total*100, 1) for b in base_order]

def get_auto_y_lim(all_pcts, step=5):
    max_val = max(all_pcts)
    y_max = math.ceil(max_val / step) * step
    return 0, y_max

def plot_combined_figure(pre_pct, mature_pct, star_pct, style):
    fig, ax = plt.subplots(figsize=style["combine_figsize"])
    bases = style["base_order"]
    n_base = len(bases)
    bw = style["bar_width"]
    bg = style["base_gap"]
    gg = style["group_gap"]

    unit = bw + bg
    x_pre = np.arange(n_base) * unit
    x_mature = x_pre + n_base * unit + gg
    x_star = x_mature + n_base * unit + gg

    for i, b in enumerate(bases):
        ax.bar(x_pre[i], pre_pct[i], bw,
               color=style["base_colors"][b],
               edgecolor=style["bar_border_color"],
               linewidth=style["bar_border_width"])
        ax.text(x_pre[i], pre_pct[i] + style["annot_offset"],
                f"{pre_pct[i]:.1f}",
                ha="center", va="bottom",
                fontsize=style["annot_fontsize"],
                fontweight=style["annot_fontweight"],
                color=style["annot_color"])

    for i, b in enumerate(bases):
        ax.bar(x_mature[i], mature_pct[i], bw,
               color=style["base_colors"][b],
               edgecolor=style["bar_border_color"],
               linewidth=style["bar_border_width"])
        ax.text(x_mature[i], mature_pct[i] + style["annot_offset"],
                f"{mature_pct[i]:.1f}",
                ha="center", va="bottom",
                fontsize=style["annot_fontsize"],
                fontweight=style["annot_fontweight"],
                color=style["annot_color"])

    for i, b in enumerate(bases):
        ax.bar(x_star[i], star_pct[i], bw,
               color=style["base_colors"][b],
               edgecolor=style["bar_border_color"],
               linewidth=style["bar_border_width"])
        ax.text(x_star[i], star_pct[i] + style["annot_offset"],
                f"{star_pct[i]:.1f}",
                ha="center", va="bottom",
                fontsize=style["annot_fontsize"],
                fontweight=style["annot_fontweight"],
                color=style["annot_color"])

    group_centers = [x_pre.mean(), x_mature.mean(), x_star.mean()]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(["Pre-miRNA", "Mature miRNA", "Star miRNA"],
                       fontsize=style["tick_fontsize"], fontweight=style["tick_fontweight"])
    ax.tick_params(axis='both', labelsize=style["tick_fontsize"], width=style["axis_line_width"])
    
    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in style["axis_spines_show"])
        ax.spines[spine].set_linewidth(style["axis_line_width"])

    if style["grid_show"]:
        ax.yaxis.grid(True, color=style["grid_color"], linewidth=style["grid_width"])
        ax.set_axisbelow(True)

    ax.set_xlabel("miRNA Type", fontsize=style["label_fontsize"], fontweight=style["label_fontweight"])
    ax.set_ylabel("Percentage (%)", fontsize=style["y_label_fontsize"], fontweight=style["label_fontweight"])
    ax.set_title("First Base Preference of miRNA Sequences",
                 fontsize=style["title_fontsize"], fontweight=style["title_fontweight"],
                 pad=style["title_pad"], loc=style["title_loc"])

    all_vals = pre_pct + mature_pct + star_pct
    y_min, y_max = get_auto_y_lim(all_vals, style["y_axis_extend_multiple"])
    ax.set_ylim(y_min, y_max)

    handles = [plt.Rectangle((0,0),1,1, fc=style["base_colors"][b],
                             ec=style["bar_border_color"],
                             lw=style["bar_border_width"]) for b in bases]
    leg = ax.legend(handles, bases,
              fontsize=style["legend_fontsize"],
              frameon=style["legend_frame_on"],
              markerscale=style["legend_markerscale"],
              loc=style["legend_loc"],
              bbox_to_anchor=style["legend_bbox_to_anchor"])
    if leg.get_frame():
        leg.get_frame().set_linewidth(style["legend_border_width"])

    plt.tight_layout()
    plt.savefig("miRNA_first_base_combined.svg", format="svg", dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def plot_single_figure(pct, seq_name, style):
    fig, ax = plt.subplots(figsize=style["single_figsize"])
    bases = style["base_order"]

    ax.bar(bases, pct,
           color=[style["base_colors"][b] for b in bases],
           edgecolor=style["bar_border_color"],
           linewidth=style["bar_border_width"])

    for i, val in enumerate(pct):
        ax.text(i, val + style["annot_offset"],
                f"{val:.1f}",
                ha="center", va="bottom",
                fontsize=style["annot_fontsize"],
                fontweight=style["annot_fontweight"],
                color=style["annot_color"])

    ax.tick_params(axis='both', labelsize=style["tick_fontsize"], width=style["axis_line_width"])
    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in style["axis_spines_show"])
        ax.spines[spine].set_linewidth(style["axis_line_width"])

    if style["grid_show"]:
        ax.yaxis.grid(True, color=style["grid_color"], linewidth=style["grid_width"])
        ax.set_axisbelow(True)

    ax.set_xlabel("First Base", fontsize=style["label_fontsize"], fontweight=style["label_fontweight"])
    ax.set_ylabel("Percentage (%)", fontsize=style["label_fontsize"], fontweight=style["label_fontweight"])
    ax.set_title(f"{seq_name} First Base Preference",
                 fontsize=style["title_fontsize"], fontweight=style["title_fontweight"],
                 pad=style["title_pad"], loc=style["title_loc"])

    y_min, y_max = get_auto_y_lim(pct, style["y_axis_extend_multiple"])
    ax.set_ylim(y_min, y_max)

    filename = f"miRNA_first_base_{seq_name.lower()}.svg"
    plt.savefig(filename, format="svg", dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def plot_alone_legend(style):
    fig, ax = plt.subplots(figsize=style["legend_figsize"])
    ax.axis("off")
    bases = style["base_order"]
    handles = [plt.Rectangle((0,0),1,1, fc=style["base_colors"][b],
                             ec=style["bar_border_color"],
                             lw=style["bar_border_width"]) for b in bases]
    leg = ax.legend(handles, bases, fontsize=style["legend_fontsize"],
                    frameon=style["legend_frame_on"], loc="center")
    if leg.get_frame():
        leg.get_frame().set_linewidth(style["legend_border_width"])
    plt.savefig("miRNA_first_base_legend.svg", format="svg", dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def main():
    args = parse_args()
    pre_bases, mature_bases, star_bases = load_first_base(args.input, CUSTOM_STYLE)
    base_order = CUSTOM_STYLE["base_order"]
    pre_pct = calc_base_percent(pre_bases, base_order)
    mature_pct = calc_base_percent(mature_bases, base_order)
    star_pct = calc_base_percent(star_bases, base_order)

    plot_combined_figure(pre_pct, mature_pct, star_pct, CUSTOM_STYLE)
    plot_single_figure(pre_pct, "Pre", CUSTOM_STYLE)
    plot_single_figure(mature_pct, "Mature", CUSTOM_STYLE)
    plot_single_figure(star_pct, "Star", CUSTOM_STYLE)
    plot_alone_legend(CUSTOM_STYLE)

    print("down!")

if __name__ == "__main__":
    main()