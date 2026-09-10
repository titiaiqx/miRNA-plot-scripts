# -*- coding: utf-8 -*-
"""
sRNA-seq 长度分布直方图
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

CUSTOM_STYLE = {
    "figure_dpi": 300,
    "font_family": "Arial",
    "figsize": (12, 5),

    "bar_color": "#4C9AC9",
    "bar_border_color": "#000000",
    "bar_border_width": 1.2,
    "bar_alpha": 1.0,
    "bar_width": 0.9,

    "annot_fontsize": 14,
    "annot_fontweight": "normal",
    "annot_color": "#000000",
    "annot_offset": 150,

    "y_axis_extend_multiple": 5,
    "axis_line_width": 1.5,
    "axis_spines_show": ["left", "bottom"],

    "grid_show": True,
    "grid_color": "#cccccc",
    "grid_width": 0.5,
    "grid_alpha": 0.5,

    "tick_fontsize": 16,
    "tick_fontweight": "normal",

    "label_fontsize": 18,
    "label_fontweight": "normal",
    "x_label": "Read Length (nt)",
    "y_label": "Total Read Count",

    "title": "sRNA-seq Total Read Length Distribution",
    "title_fontsize": 16,
    "title_fontweight": "bold",
    "title_pad": 20,
    "title_loc": "center",

    "save_formats": ["svg", "png"],
}
# ======================================================================================

plt.rcParams['font.family'] = CUSTOM_STYLE["font_family"]
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

def parse_args():
    parser = argparse.ArgumentParser(description="sRNA-seq length distribution histogram")
    parser.add_argument("-i", "--input", required=True, help="输入统计文件")
    parser.add_argument("-o", "--output", default="total_length_dist", help="输出文件名")
    return parser.parse_args()

def load_length_data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df[df["Length"] != "TotalReads"].copy()
    df["Length"] = df["Length"].astype(int)
    df["Total"] = df["Total"].astype(int)
    df = df.sort_values("Length")
    return df

def get_auto_y_lim(values, step=5):
    max_val = max(values)
    y_max = math.ceil(max_val / step) * step
    return 0, y_max

def plot_length_distribution(df, style, output_prefix):
    fig, ax = plt.subplots(figsize=style["figsize"])
    
    x = df["Length"].values
    y = df["Total"].values

    ax.bar(
        x, y,
        width=style["bar_width"],
        color=style["bar_color"],
        edgecolor=style["bar_border_color"],
        linewidth=style["bar_border_width"],
        alpha=style["bar_alpha"]
    )

    # for xi, yi in zip(x, y):
    #     ax.text(
    #         xi, yi + style["annot_offset"],
    #         f"{yi}",
    #         ha="center", va="bottom",
    #         fontsize=style["annot_fontsize"],
    #         fontweight=style["annot_fontweight"],
    #         color=style["annot_color"]
    #     )

    ax.set_xlabel(
        style["x_label"],
        fontsize=style["label_fontsize"],
        fontweight=style["label_fontweight"]
    )
    ax.set_ylabel(
        style["y_label"],
        fontsize=style["label_fontsize"],
        fontweight=style["label_fontweight"]
    )

    ax.tick_params(
        axis='both',
        labelsize=style["tick_fontsize"],
        width=style["axis_line_width"]
    )

    for spine in ax.spines:
        ax.spines[spine].set_visible(spine in style["axis_spines_show"])
        ax.spines[spine].set_linewidth(style["axis_line_width"])
        ax.spines[spine].set_color(style["bar_border_color"])

    if style["grid_show"]:
        ax.yaxis.grid(
            True,
            color=style["grid_color"],
            linewidth=style["grid_width"],
            alpha=style["grid_alpha"]
        )
        ax.set_axisbelow(True)

    ax.set_title(
        style["title"],
        fontsize=style["title_fontsize"],
        fontweight=style["title_fontweight"],
        pad=style["title_pad"],
        loc=style["title_loc"]
    )

    y_min, y_max = get_auto_y_lim(y, step=style["y_axis_extend_multiple"])
    ax.set_ylim(y_min, y_max)

    plt.tight_layout()

    for fmt in style["save_formats"]:
        plt.savefig(
            f"{output_prefix}.{fmt}",
            format=fmt,
            dpi=style["figure_dpi"],
            bbox_inches="tight"
        )
    plt.close()

def main():
    args = parse_args()
    df = load_length_data(args.input)
    plot_length_distribution(df, CUSTOM_STYLE, args.output)
    print(f"down!")

if __name__ == "__main__":
    main()