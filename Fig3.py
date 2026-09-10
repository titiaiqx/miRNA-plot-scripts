# -*- coding: utf-8 -*-
"""
差异分析火山图
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

POINT_SIZE = 13
ALPHA = 0.8
UP_COLOR = "#E2745E"
DOWN_COLOR = "#4C9AC9"
NS_COLOR = "#B2B2B2"

EDGE_WIDTH_UP = 0.5
EDGE_COLOR_UP = "#9A1D05"

EDGE_WIDTH_DOWN = 0.5
EDGE_COLOR_DOWN = "#004488"

EDGE_WIDTH_NS = 0.5
EDGE_COLOR_NS = "#434343"

THRESHOLD_LINE_WIDTH = 1.0
THRESHOLD_COLOR = "#5d5d5d"
DASH_PATTERN = (4, 2)

SHOW_THRESHOLD_LABELS = False
LABEL_FONTSIZE = 9
LABEL_COLOR = "#111111"

AXIS_WIDTH = 2.0

fc_threshold = 2
p_threshold = 0.05

def plot_volcano(input_file, use_extended_y=False):
    df = pd.read_csv(input_file, sep='\t')
    
    x = df['log2FoldChange']
    pvals = df['padj']
    pvals = np.clip(pvals, 1e-300, 1.0)
    y = -np.log10(pvals)
    
    up = (x >= fc_threshold) & (pvals <= p_threshold)
    down = (x <= -fc_threshold) & (pvals <= p_threshold)
    ns = ~(up | down)

    plt.figure(figsize=(6, 6), dpi=300)
    ax = plt.gca()

    ax.scatter(x[ns], y[ns],
               c=NS_COLOR, s=POINT_SIZE, alpha=ALPHA,
               edgecolors=EDGE_COLOR_NS, linewidths=EDGE_WIDTH_NS)
    ax.scatter(x[down], y[down],
               c=DOWN_COLOR, s=POINT_SIZE, alpha=ALPHA,
               edgecolors=EDGE_COLOR_DOWN, linewidths=EDGE_WIDTH_DOWN)
    ax.scatter(x[up], y[up],
               c=UP_COLOR, s=POINT_SIZE, alpha=ALPHA,
               edgecolors=EDGE_COLOR_UP, linewidths=EDGE_WIDTH_UP)

    threshold_y = -np.log10(p_threshold)
    ax.axvline(x=fc_threshold,
               color=THRESHOLD_COLOR, linestyle='--',
               linewidth=THRESHOLD_LINE_WIDTH, dashes=DASH_PATTERN)
    ax.axvline(x=-fc_threshold,
               color=THRESHOLD_COLOR, linestyle='--',
               linewidth=THRESHOLD_LINE_WIDTH, dashes=DASH_PATTERN)
    ax.axhline(y=threshold_y,
               color=THRESHOLD_COLOR, linestyle='--',
               linewidth=THRESHOLD_LINE_WIDTH, dashes=DASH_PATTERN)

    if SHOW_THRESHOLD_LABELS:
        ax.text(fc_threshold + 0.09, ax.get_ylim()[1] * 0.94,
                f"{fc_threshold}", fontsize=LABEL_FONTSIZE, color=LABEL_COLOR)
        ax.text(-fc_threshold - 0.5, ax.get_ylim()[1] * 0.94,
                f"-{fc_threshold}", fontsize=LABEL_FONTSIZE, color=LABEL_COLOR)
        ax.text(ax.get_xlim()[1] * 0.81, threshold_y + 0.25,
                f"{threshold_y:.2f}", fontsize=LABEL_FONTSIZE, color=LABEL_COLOR)

    if use_extended_y:
        ax.set_yscale('function', functions=(
            lambda val: np.log10(val + 1),
            lambda val: 10**val - 1
        ))
        ax.set_ylim(bottom=0)
    else:
        ax.set_ylim(bottom=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(AXIS_WIDTH)
    ax.spines['bottom'].set_linewidth(AXIS_WIDTH)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend('', frameon=False)

    plt.tight_layout(pad=0.1)
    plt.savefig("volcano_plot.png", dpi=300, bbox_inches='tight')
    plt.savefig("volcano_plot.svg", bbox_inches='tight')
    plt.close()

    print("down!")
    if use_extended_y:
        print("open")
    if SHOW_THRESHOLD_LABELS:
        print("show")
    print(f"up: {sum(up)}")
    print(f"down: {sum(down)}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='火山图绘制')
    parser.add_argument('-i', '--input', required=True, help='tsv文件')
    parser.add_argument('-ex', '--extend', action='store_true', help='纵轴下密上疏')
    args = parser.parse_args()
    plot_volcano(args.input, use_extended_y=args.extend)