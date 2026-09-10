# -*- coding: utf-8 -*-
"""
miRNA upset图
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

def load_and_preprocess(input_path):

    df = pd.read_csv(input_path)

    groups = {
        "FL": ["FL-1", "FL-2", "FL-3"],
        "FS": ["FS-1", "FS-2", "FS-3"],
        "ML": ["ML-1", "ML-2", "ML-3"],
        "MS": ["MS-1", "MS-2", "MS-3"]
    }

    for group, cols in groups.items():
        df[group] = df[cols].sum(axis=1)

    df["total_count"] = df[list(groups.keys())].sum(axis=1)
    df = df.sort_values(by="total_count", ascending=False).reset_index(drop=True)

    indicator_df = pd.DataFrame()
    for group in groups.keys():
        indicator_df[group] = df[group] > 0

    indicator_df["gene_id"] = df["gene_id"]
    indicator_df = indicator_df.set_index("gene_id")
    return indicator_df

def set_custom_style():

    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9

    COLORS = {
        "facecolor": "#4C9AC9",
        "other_dots_color": "#bcc9d1",
        "shading_color": "#EBEBEB",
        "edgecolor": "#000000",
        "totals_color": "#E2745E"
    }

    LINE_PARAMS = {
        "linewidth": 0.8,
        "linestyle": "-",
        "with_lines": True
    }
    LAYOUT_PARAMS = {
        "figsize": (10, 6),
        "element_size": 25,
        "intersection_plot_elements": 8,
        "totals_plot_elements": 4,
        "sort_by": "cardinality",
        "sort_categories_by": "cardinality"
    }

    LABEL_PARAMS = {
        "show_counts": True,
        "count_format": "{:d}",
        "show_percentages": False
    }

    return COLORS, LINE_PARAMS, LAYOUT_PARAMS, LABEL_PARAMS

def plot_upset(indicator_df):

    COLORS, LINE_PARAMS, LAYOUT_PARAMS, LABEL_PARAMS = set_custom_style()
    upset_data = from_indicators(indicator_df.columns, data=indicator_df)
    fig = plt.figure(figsize=LAYOUT_PARAMS["figsize"])

    upset = UpSet(
        upset_data,
        orientation="horizontal",
        sort_by=LAYOUT_PARAMS["sort_by"],
        sort_categories_by=LAYOUT_PARAMS["sort_categories_by"],
        facecolor=COLORS["facecolor"],
        other_dots_color=COLORS["other_dots_color"],
        shading_color=COLORS["shading_color"],
        with_lines=LINE_PARAMS["with_lines"],
        element_size=LAYOUT_PARAMS["element_size"],
        intersection_plot_elements=LAYOUT_PARAMS["intersection_plot_elements"],
        totals_plot_elements=LAYOUT_PARAMS["totals_plot_elements"],
        show_counts=LABEL_PARAMS["show_counts"],
        show_percentages=LABEL_PARAMS["show_percentages"]
    )
    subplots = upset.plot(fig=fig)

    ax_intersections = subplots["intersections"]
    ax_intersections.set_ylim(0, 50)   # 截断
    ax_intersections.set_yticks([0,10,20,30,40,50])

    ax_totals = subplots["totals"]
    for bar in ax_totals.patches:
        bar.set_facecolor(COLORS["totals_color"])
        bar.set_edgecolor(COLORS["edgecolor"])
        bar.set_linewidth(LINE_PARAMS["linewidth"])

    for bar in ax_intersections.patches:
        bar.set_edgecolor(COLORS["edgecolor"])
        bar.set_linewidth(LINE_PARAMS["linewidth"])

    ax_matrix = subplots["matrix"]
    for line in ax_matrix.lines:
        line.set_color(COLORS["edgecolor"])
        line.set_linewidth(LINE_PARAMS["linewidth"])
        line.set_linestyle(LINE_PARAMS["linestyle"])

    plt.suptitle(
        "miRNA Presence in Four Sample Groups (FL, FS, ML, MS)",
        fontsize=plt.rcParams['axes.titlesize'],
        y=0.98
    )
    plt.tight_layout()

    plt.savefig(
        "miRNA_upset_plot.svg",
        bbox_inches="tight",
        facecolor="white"
    )
    plt.savefig(
        "miRNA_upset_plot.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white"
    )
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Generate Upset Plot for miRNA Expression")
    parser.add_argument("-i", "--input", required=True, help="Input miRNA count CSV file path")
    args = parser.parse_args()
    indicator_df = load_and_preprocess(args.input)
    plot_upset(indicator_df)
    print("down!\n1. miRNA_upset_plot.svg\n2. miRNA_upset_plot.png")

if __name__ == "__main__":
    main()