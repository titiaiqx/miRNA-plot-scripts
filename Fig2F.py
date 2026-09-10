# -*- coding: utf-8 -*-
"""
miRNA表达箱线图 + 散点图 + 直方图
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

BOX_PLOT_SETTINGS = {
    "fig_width": 9,
    "fig_height": 7,
    "dpi": 300,
    "font": "Arial",

    "title_text": "",
    "title_size": 26,
    "title_pad": 15,

    "y_label_text": "Log10(TPM)",
    "y_label_size": 26,
    "x_label_text": "",

    "tick_label_size": 22,
    "tick_width": 2,

    "axis_line_width": 2,
    "show_left_spine": True,
    "show_bottom_spine": True,
    "show_right_spine": False,
    "show_top_spine": False,

    "show_grid": True,
    "grid_color": "#EAEAEA",
    "grid_width": 1.2,

    "y_min": -2,
    "y_max": 5,

    "box_line_width": 2,
    "box_width": 0.7,
    "box_alpha": 0.7,

    "scatter_size": 3,
    "scatter_alpha": 0.8,
    "scatter_jitter": 0.22,
    "max_scatter_points": 4000,

    "save_path": "miRNA_boxplot_strip.svg"
}

HIST_PLOT_SETTINGS = {
    "fig_width": 10,
    "fig_height": 6,
    "dpi": 300,
    "font": "Arial",

    "title_text": "miRNA Expression Distribution",
    "title_size": 24,
    "title_pad": 15,

    "x_label_text": "Log10(TPM)",
    "y_label_text": "Count",
    "label_size": 22,

    "tick_label_size": 20,
    "tick_width": 2,

    "axis_line_width": 2,
    "show_left_spine": True,
    "show_bottom_spine": True,
    "show_right_spine": False,
    "show_top_spine": False,

    "show_grid": True,
    "grid_color": "#EAEAEA",
    "grid_width": 1.2,

    "bin_count": 60,
    "bar_border_color": "black",
    "bar_border_width": 0,
    "bar_alpha": 0.9,

    "color_start": "#2E86AB",
    "color_end": "#2E86AB",

    "save_path": "miRNA_expression_histogram.svg"
}

GROUP_COLORS = {
    "FL": "#72B4DE",
    "FS": "#c6ddeb",
    "ML": "#F49A88",
    "MS": "#FDD2BE"
}

BOX_EDGE_COLORS = {
    "FL": "#000000",
    "FS": "#000000",
    "ML": "#000000",
    "MS": "#000000"
}

SCATTER_DARK_COLORS = {
    "FL": "#07659F",
    "FS": "#6194b4",
    "ML": "#B03820",
    "MS": "#D4774B"
}

GROUP_ORDER = ["FL", "FS", "ML", "MS"]

def parse_args():
    parser = argparse.ArgumentParser(description="箱线散点图 + 直方图")
    parser.add_argument("-i", required=True, help="miRNA TPM 表达矩阵 TSV 文件")
    return parser.parse_args()

def load_data(file_path):
    return pd.read_csv(file_path, sep="\t", index_col=0)

def get_sample_groups(sample_list):
    groups = {}
    for s in sample_list:
        if s.startswith("FL-"): groups[s] = "FL"
        elif s.startswith("FS-"): groups[s] = "FS"
        elif s.startswith("ML-"): groups[s] = "ML"
        elif s.startswith("MS-"): groups[s] = "MS"
    return groups

def plot_box_scatter(df, groups, cfg):
    melt = df.melt(var_name="Sample", value_name="TPM")
    melt["Group"] = melt["Sample"].map(groups)
    melt = melt[melt["TPM"] > 0].copy()
    melt["log10_TPM"] = np.log10(melt["TPM"])

    plt.rcParams['font.family'] = cfg["font"]
    fig, ax = plt.subplots(figsize=(cfg["fig_width"], cfg["fig_height"]), dpi=cfg["dpi"])

    group_data_list = []
    for g in GROUP_ORDER:
        sub = melt[melt["Group"] == g]["log10_TPM"].values
        group_data_list.append(sub)

    bp = ax.boxplot(
        group_data_list,
        positions=list(range(len(GROUP_ORDER))),
        patch_artist=True,
        widths=cfg["box_width"],
        showfliers=False,
        labels=GROUP_ORDER
    )

    for idx, g_name in enumerate(GROUP_ORDER):
        fill_c = GROUP_COLORS[g_name]
        edge_c = BOX_EDGE_COLORS[g_name]
        box = bp["boxes"][idx]
        box.set_facecolor(fill_c)
        box.set_edgecolor(edge_c)
        box.set_linewidth(cfg["box_line_width"])
        box.set_alpha(cfg["box_alpha"])
        med = bp["medians"][idx]
        med.set_color(edge_c)
        med.set_linewidth(cfg["box_line_width"])
        whisker1 = bp["whiskers"][idx * 2]
        whisker2 = bp["whiskers"][idx * 2 + 1]
        whisker1.set_color(edge_c)
        whisker2.set_color(edge_c)
        whisker1.set_linewidth(cfg["box_line_width"])
        whisker2.set_linewidth(cfg["box_line_width"])
        cap1 = bp["caps"][idx * 2]
        cap2 = bp["caps"][idx * 2 + 1]
        cap1.set_color(edge_c)
        cap2.set_color(edge_c)
        cap1.set_linewidth(cfg["box_line_width"])
        cap2.set_linewidth(cfg["box_line_width"])

    sub = melt.sample(n=min(cfg["max_scatter_points"], len(melt)), random_state=42)
    sns.stripplot(
        x="Group", y="log10_TPM", hue="Group",
        data=sub, order=GROUP_ORDER,
        palette=SCATTER_DARK_COLORS, ax=ax, legend=False,
        size=cfg["scatter_size"],
        alpha=cfg["scatter_alpha"],
        jitter=cfg["scatter_jitter"]
    )


    ax.set_title(cfg["title_text"], fontsize=cfg["title_size"], pad=cfg["title_pad"])
    ax.set_ylabel(cfg["y_label_text"], fontsize=cfg["y_label_size"])
    ax.set_xlabel(cfg["x_label_text"])
    ax.tick_params(axis="both", labelsize=cfg["tick_label_size"], width=cfg["tick_width"])

    ax.spines["left"].set_visible(cfg["show_left_spine"])
    ax.spines["bottom"].set_visible(cfg["show_bottom_spine"])
    ax.spines["right"].set_visible(cfg["show_right_spine"])
    ax.spines["top"].set_visible(cfg["show_top_spine"])
    for s in ax.spines.values():
        s.set_linewidth(cfg["axis_line_width"])

    if cfg["show_grid"]:
        ax.yaxis.grid(True, color=cfg["grid_color"], linewidth=cfg["grid_width"])
        ax.set_axisbelow(True)

    ax.set_ylim(cfg["y_min"], cfg["y_max"])

    plt.tight_layout()
    plt.savefig(cfg["save_path"], bbox_inches="tight")
    plt.close()

def plot_histogram(df, groups, cfg):
    melt = df.melt(var_name="Sample", value_name="TPM")
    melt["Group"] = melt["Sample"].map(groups)
    melt = melt[melt["TPM"] > 0].copy()
    melt["log10_TPM"] = np.log10(melt["TPM"])

    plt.rcParams['font.family'] = cfg["font"]
    fig, ax = plt.subplots(figsize=(cfg["fig_width"], cfg["fig_height"]), dpi=cfg["dpi"])

    cmap = LinearSegmentedColormap.from_list("custom", [cfg["color_start"], cfg["color_end"]], N=256)
    n, bins, patches = ax.hist(
        melt["log10_TPM"], bins=cfg["bin_count"],
        edgecolor=cfg["bar_border_color"],
        linewidth=cfg["bar_border_width"]
    )

    norm = plt.Normalize(melt["log10_TPM"].min(), melt["log10_TPM"].max())
    for p, x in zip(patches, np.linspace(melt["log10_TPM"].min(), melt["log10_TPM"].max(), len(patches))):
        p.set_facecolor(cmap(norm(x)))
        p.set_alpha(cfg["bar_alpha"])

    ax.set_title(cfg["title_text"], fontsize=cfg["title_size"], pad=cfg["title_pad"])
    ax.set_xlabel(cfg["x_label_text"], fontsize=cfg["label_size"])
    ax.set_ylabel(cfg["y_label_text"], fontsize=cfg["label_size"])
    ax.tick_params(axis="both", labelsize=cfg["tick_label_size"], width=cfg["tick_width"])

    ax.spines["left"].set_visible(cfg["show_left_spine"])
    ax.spines["bottom"].set_visible(cfg["show_bottom_spine"])
    ax.spines["right"].set_visible(cfg["show_right_spine"])
    ax.spines["top"].set_visible(cfg["show_top_spine"])
    for s in ax.spines.values():
        s.set_linewidth(cfg["axis_line_width"])

    if cfg["show_grid"]:
        ax.yaxis.grid(True, color=cfg["grid_color"], linewidth=cfg["grid_width"])
        ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(cfg["save_path"], bbox_inches="tight")
    plt.close()

def main():
    args = parse_args()
    df = load_data(args.i)
    sample_groups = get_sample_groups(df.columns)

    plot_box_scatter(df, sample_groups, BOX_PLOT_SETTINGS)
    plot_histogram(df, sample_groups, HIST_PLOT_SETTINGS)
    print("down!")

if __name__ == "__main__":
    main()