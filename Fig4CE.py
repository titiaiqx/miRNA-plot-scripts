# -*- coding: utf-8 -*-
"""
miRNA-靶基因散点图
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

CUSTOM_STYLE = {

    "figure_dpi": 300,
    "font_family": "Arial",
    "figsize": (8, 8),
    "legend_figsize": (5, 4),

    "base_color": "#D3D3D3",
    "base_alpha": 0.4,
    "base_size": 100,

    "q2_color": "#457B9D",
    "q4_color": "#E63946",
    "highlight_size": 150,
    "edge_color": "black",
    "edge_width": 1.5,
    "pl_alpha": 1.0,

    "mirna_threshold": 1.0,
    "gene_threshold": 1.0,

    "zero_vline_color": "#525252",
    "zero_hline_color": "#525252",
    "zero_line_width": 1.4,
    "zero_line_style": "solid",

    "thresh_line_color": "#999999",
    "thresh_line_width": 1.2,
    "thresh_line_style": "--",

    "axis_line_width": 2,
    "axis_spines_show": ["left", "bottom", "top", "right"],
    "x_lim": (-15, 15),
    "y_lim": (-15, 15),

    "grid_show": True,
    "grid_color": "#EAEAEA",
    "grid_width": 1.2,
    "grid_alpha": 0.8,

    "tick_fontsize": 18,
    "tick_fontweight": "normal",

    "label_fontsize": 20,
    "label_fontweight": "normal",
    "x_label": "miRNA log2FC",
    "y_label": "Target gene log2FC",

    "title_show": False,
    "title_text": "miRNA-Target Gene Expression Quadrant Plot",
    "title_fontsize": 18,
    "title_fontweight": "bold",
    "title_pad": 20,

    "legend_fontsize": 12,
    "legend_fontweight": "normal",
    "legend_border_width": 1.2,
    "legend_frame_on": True,
    "legend_markerscale": 1,

    "save_format": "svg",
    "main_plot_name": "mirna_target_quadrant_plot",
    "legend_plot_name": "mirna_target_legend_only",
    "final_filtered_output": "final_q2_q4_significant_pairs.txt",
}

plt.rcParams['font.family'] = CUSTOM_STYLE["font_family"]
plt.rcParams['axes.unicode_minus'] = False


def parse_args():
    parser = argparse.ArgumentParser(description="miRNA-靶基因散点图")
    parser.add_argument("-p", required=True, help="miRNA-靶基因")
    parser.add_argument("-pl", required=True, help="显著配对")
    return parser.parse_args()


def load_data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.dropna(subset=["miRNA_log2FC", "gene_log2FC"]).copy()
    df["miRNA_log2FC"] = df["miRNA_log2FC"].astype(float)
    df["gene_log2FC"] = df["gene_log2FC"].astype(float)
    return df


def save_final_q2q4_pairs(df_pl, mirna_cut, gene_cut, out_path):

    q2_mask = (df_pl["miRNA_log2FC"] > mirna_cut) & (df_pl["gene_log2FC"] < -gene_cut)
    q4_mask = (df_pl["miRNA_log2FC"] < -mirna_cut) & (df_pl["gene_log2FC"] > gene_cut)
    total_mask = q2_mask | q4_mask

    final_df = df_pl[total_mask].copy()
    final_df.to_csv(out_path, sep="\t", index=False)
    print(f"all: {len(final_df)}")
    print(f"2: {sum(q2_mask[total_mask])}")
    print(f"4: {sum(q4_mask[total_mask])}")


def plot_main_figure(style):
    fig, ax = plt.subplots(figsize=style["figsize"])
    args = parse_args()
    df_all = load_data(args.p)
    df_pl = load_data(args.pl)

    mirna_cut = style["mirna_threshold"]
    gene_cut = style["gene_threshold"]

    save_final_q2q4_pairs(df_pl, mirna_cut, gene_cut, style["final_filtered_output"])

    x_all = df_all["miRNA_log2FC"]
    y_all = df_all["gene_log2FC"]

    ax.scatter(x_all, y_all, c=style["base_color"], alpha=style["base_alpha"],
               s=style["base_size"], zorder=1)

    q2_all = (x_all > mirna_cut) & (y_all < -gene_cut)
    q4_all = (x_all < -mirna_cut) & (y_all > gene_cut)
    ax.scatter(x_all[q2_all], y_all[q2_all], c=style["q2_color"],
               alpha=style["base_alpha"], s=style["highlight_size"], zorder=2)
    ax.scatter(x_all[q4_all], y_all[q4_all], c=style["q4_color"],
               alpha=style["base_alpha"], s=style["highlight_size"], zorder=2)

    x_pl = df_pl["miRNA_log2FC"]
    y_pl = df_pl["gene_log2FC"]
    pl_q2 = (x_pl > mirna_cut) & (y_pl < -gene_cut)
    pl_q4 = (x_pl < -mirna_cut) & (y_pl > gene_cut)

    ax.scatter(x_pl[pl_q2], y_pl[pl_q2], c=style["q2_color"],
               edgecolors=style["edge_color"], linewidth=style["edge_width"],
               alpha=style["pl_alpha"], s=style["highlight_size"], zorder=5)
    ax.scatter(x_pl[pl_q4], y_pl[pl_q4], c=style["q4_color"],
               edgecolors=style["edge_color"], linewidth=style["edge_width"],
               alpha=style["pl_alpha"], s=style["highlight_size"], zorder=5)

    ax.axvline(0, color=style["zero_vline_color"], linestyle=style["zero_line_style"],
               linewidth=style["zero_line_width"], zorder=100)
    ax.axhline(0, color=style["zero_hline_color"], linestyle=style["zero_line_style"],
               linewidth=style["zero_line_width"], zorder=100)

    ax.axvline(x=mirna_cut, color=style["thresh_line_color"],
               linestyle=style["thresh_line_style"], linewidth=style["thresh_line_width"])
    ax.axvline(x=-mirna_cut, color=style["thresh_line_color"],
               linestyle=style["thresh_line_style"], linewidth=style["thresh_line_width"])
    ax.axhline(y=gene_cut, color=style["thresh_line_color"],
               linestyle=style["thresh_line_style"], linewidth=style["thresh_line_width"])
    ax.axhline(y=-gene_cut, color=style["thresh_line_color"],
               linestyle=style["thresh_line_style"], linewidth=style["thresh_line_width"])
    
    ax.set_xlim(style["x_lim"])
    ax.set_ylim(style["y_lim"])

    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in style["axis_spines_show"])
        ax.spines[spine].set_linewidth(style["axis_line_width"])

    if style["grid_show"]:
        ax.grid(True, color=style["grid_color"], linewidth=style["grid_width"], alpha=style["grid_alpha"])
        ax.set_axisbelow(True)

    ax.tick_params(axis="both", labelsize=style["tick_fontsize"], width=style["axis_line_width"])
    ax.set_xlabel(style["x_label"], fontsize=style["label_fontsize"], fontweight=style["label_fontweight"])
    ax.set_ylabel(style["y_label"], fontsize=style["label_fontsize"], fontweight=style["label_fontweight"])

    if style["title_show"]:
        ax.set_title(style["title_text"], fontsize=style["title_fontsize"],
                     fontweight=style["title_fontweight"], pad=style["title_pad"])

    plt.tight_layout()
    main_out = f"{style['main_plot_name']}.{style['save_format']}"
    plt.savefig(main_out, dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def plot_only_legend(style):
    fig, ax = plt.subplots(figsize=style["legend_figsize"])
    ax.axis("off")

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=style["base_color"],
               markersize=10, alpha=style["base_alpha"], label='All miRNA-target pairs'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=style["q2_color"],
               markersize=10, label='miRNA up / Target down (Q2)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=style["q4_color"],
               markersize=10, label='miRNA down / Target up (Q4)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=style["q2_color"],
               markeredgecolor=style["edge_color"], markeredgewidth=style["edge_width"],
               markersize=10, label='Significant Q2/Q4 pairs'),
    ]

    leg = ax.legend(handles=legend_elements, loc="center",
                    fontsize=style["legend_fontsize"],
                    frameon=style["legend_frame_on"],
                    markerscale=style["legend_markerscale"])

    if leg.get_frame():
        leg.get_frame().set_linewidth(style["legend_border_width"])

    plt.tight_layout()
    leg_out = f"{style['legend_plot_name']}.{style['save_format']}"
    plt.savefig(leg_out, dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def main():
    plot_main_figure(CUSTOM_STYLE)
    plot_only_legend(CUSTOM_STYLE)
    print("down!")


if __name__ == "__main__":
    main()