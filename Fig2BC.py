#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
miRNA长度分布
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter

CUSTOM_STYLE = {
    "figure_dpi": 300,
    "font_family": "Arial",

    "pre_fig_size": (13, 7),
    "pre_group_step": 50,
    "pre_max_length": 300,
    "pre_gradient_start": "#2E86AB",
    "pre_gradient_end": "#A23B72",
    "pre_border_color": "black",
    "pre_border_width": 2,
    "pre_curve_color": "#54BD12",
    "pre_curve_width": 5,
    "pre_tick_fontsize": 16,
    "pre_label_fontsize": 18,
    "pre_title_fontsize": 18,
    "pre_axis_line_width": 2,
    "pre_grid_show": True,
    "pre_grid_color": "#EAEAEA",
    "pre_grid_width": 1.5,
    "y_tick_step": 50,

    "mature_fig_size": (10, 6),
    "mature_bar_width": 0.6,
    "mature_gradient_start": "#2E86AB",
    "mature_gradient_end": "#A23B72",
    "mature_border_color": "black",
    "mature_border_width": 2,
    "mature_tick_fontsize": 16,
    "mature_label_fontsize": 18,
    "mature_title_fontsize": 18,
    "mature_axis_line_width": 2,
    "mature_grid_show": True,
    "mature_grid_color": "#EAEAEA",
    "mature_grid_width": 1.5,

    "annot_fontsize": 17,
    "pre_col": 11,
    "mature_col": 18,
    "star_col": 20,
}

plt.rcParams['font.family'] = CUSTOM_STYLE["font_family"]
plt.rcParams['axes.unicode_minus'] = False

def parse_args():
    parser = argparse.ArgumentParser(description="miRNA长度分布分析")
    parser.add_argument("-i", "--input", required=True, help="tsv文件")
    return parser.parse_args()

def load_data(file_path, style):
    df = pd.read_csv(file_path, sep="\t", header=None)
    pre_seq = df.iloc[:, style["pre_col"]].dropna().astype(str).str.strip()
    mature_seq = df.iloc[:, style["mature_col"]].dropna().astype(str).str.strip()
    star_seq = df.iloc[:, style["star_col"]].dropna().astype(str).str.strip()
    pre_seq = [s for s in pre_seq if s]
    mature_seq = [s for s in mature_seq if s]
    star_seq = [s for s in star_seq if s]
    print(f"Pre-miRNA: {len(pre_seq)}")
    print(f"Mature: {len(mature_seq)}")
    print(f"Star: {len(star_seq)}")
    return pre_seq, mature_seq, star_seq

def get_pre_hist_data(sequences, style):
    lengths = [len(s) for s in sequences]
    bins = np.arange(0, style["pre_max_length"] + style["pre_group_step"], style["pre_group_step"])
    counts, edges = np.histogram(lengths, bins=bins)
    return edges, counts

def get_mature_hist_data(sequences):
    lengths = [len(s) for s in sequences]
    cnt = sorted(Counter(lengths).items())
    return ([x for x, y in cnt], [y for x, y in cnt]) if cnt else ([], [])

def get_cumulative_curve(sequences):
    lengths = [len(s) for s in sequences]
    if not lengths:
        return np.array([0]), np.array([0])
    x_raw = np.sort(lengths)
    y_raw = np.arange(1, len(x_raw)+1) / len(x_raw) * 100
    x_smooth = np.linspace(x_raw.min(), x_raw.max(), 2000)
    try:
        spl = make_interp_spline(x_raw, y_raw, k=2, bc_type="natural")
        return x_smooth, spl(x_smooth)
    except:
        return x_raw, y_raw

def get_regular_y_ticks(max_val, step):
    end = int(np.ceil(max_val / step) * step)
    return np.arange(0, end + step, step)

def plot_pre_hist_with_curve(bin_edges, counts, cum_x, cum_y, style):
    fig, ax1 = plt.subplots(figsize=style["pre_fig_size"])
    n = len(counts)
    cmap = LinearSegmentedColormap.from_list("pre", [style["pre_gradient_start"], style["pre_gradient_end"]])
    colors = cmap(np.linspace(0, 1, n))

    bars = ax1.bar(
        bin_edges[:-1], counts, width=style["pre_group_step"],
        color=colors, edgecolor=style["pre_border_color"],
        linewidth=style["pre_border_width"], align="edge"
    )

    for bar, c in zip(bars, counts):
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width()/2., height + 0.1,
            f'{int(c)}', ha='center', va='bottom',
            fontsize=style["annot_fontsize"]
        )

    ax2 = ax1.twinx()

    ax2.spines['right'].set_linewidth(style["pre_axis_line_width"])
    ax2.spines['top'].set_visible(False)
    ax2.tick_params(axis='y', width=style["pre_axis_line_width"])

    ax2.plot(cum_x, cum_y, color=style["pre_curve_color"], linewidth=style["pre_curve_width"])
    ax2.set_ylim(0, 105)

    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_linewidth(style["pre_axis_line_width"])
    ax1.spines['bottom'].set_linewidth(style["pre_axis_line_width"])

    if style["pre_grid_show"]:
        ax1.yaxis.grid(True, color=style["pre_grid_color"], linewidth=style["pre_grid_width"])
        ax1.set_axisbelow(True)

    y_max = ax1.get_ylim()[1]
    y_ticks = get_regular_y_ticks(y_max, style["y_tick_step"])
    ax1.set_yticks(y_ticks)
    ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{int(v)}"))

    ax1.tick_params(axis='y', labelsize=style["pre_tick_fontsize"])
    ax2.tick_params(axis='y', labelsize=style["pre_tick_fontsize"])
    ax1.tick_params(axis='x', labelsize=style["pre_tick_fontsize"])

    ax1.set_xlabel("Length (nt)", fontsize=style["pre_label_fontsize"])
    ax1.set_ylabel("Count", fontsize=style["pre_label_fontsize"])
    ax2.set_ylabel("Cumulative %", fontsize=style["pre_label_fontsize"])
    plt.title("Pre-miRNA Length Distribution", fontsize=style["pre_title_fontsize"], pad=20)

    plt.savefig("miRNA_pre_length_dist.svg", format="svg", dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def plot_mature_hist(xs, ys, style):
    if not xs:
        return
    fig, ax = plt.subplots(figsize=style["mature_fig_size"])
    n = len(xs)
    cmap = LinearSegmentedColormap.from_list("mat", [style["mature_gradient_start"], style["mature_gradient_end"]])
    colors = cmap(np.linspace(0, 1, n))

    bars = ax.bar(
        xs, ys, width=style["mature_bar_width"],
        color=colors, edgecolor=style["mature_border_color"],
        linewidth=style["mature_border_width"]
    )

    for bar, y in zip(bars, ys):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2., height + 0.1,
            f'{int(y)}', ha='center', va='bottom',
            fontsize=style["annot_fontsize"]
        )

    ax.set_xticks(np.arange(int(min(xs)), int(max(xs)) + 1, 1))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v)}"))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(style["mature_axis_line_width"])
    ax.spines['bottom'].set_linewidth(style["mature_axis_line_width"])

    if style["mature_grid_show"]:
        ax.yaxis.grid(True, color=style["mature_grid_color"], linewidth=style["mature_grid_width"])
        ax.set_axisbelow(True)

    y_max = ax.get_ylim()[1]
    y_ticks = get_regular_y_ticks(y_max, style["y_tick_step"])
    ax.set_yticks(y_ticks)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{int(v)}"))

    ax.tick_params(axis='y', labelsize=style["mature_tick_fontsize"])
    ax.tick_params(axis='x', labelsize=style["mature_tick_fontsize"])

    ax.set_xlabel("Length (nt)", fontsize=style["mature_label_fontsize"])
    ax.set_ylabel("Count", fontsize=style["mature_label_fontsize"])
    plt.title("Mature miRNA Length Distribution", fontsize=style["mature_title_fontsize"], pad=20)

    plt.savefig("miRNA_mature_length_dist.svg", format="svg", dpi=style["figure_dpi"], bbox_inches="tight")
    plt.close()

def save_statistics(pre_seq, mature_seq, star_seq, bin_edges, pre_counts, mature_x, mature_y):
    stats = []
    stats.append(["Type", "Total"])
    stats.append(["Pre-miRNA", len(pre_seq)])
    stats.append(["Mature", len(mature_seq)])
    stats.append(["Star", len(star_seq)])
    stats.append(["", ""])

    def len_stats(name, seqs):
        if not seqs: return
        ls = [len(s) for s in seqs]
        stats.append([f"{name} Length", ""])
        stats.append(["Min", min(ls)])
        stats.append(["Max", max(ls)])
        stats.append(["Mean", round(np.mean(ls), 2)])
        stats.append(["Median", round(np.median(ls), 2)])
        stats.append(["", ""])

    len_stats("Pre", pre_seq)
    len_stats("Mature", mature_seq)
    len_stats("Star", star_seq)

    stats.append(["Pre Range", "Count"])
    for i in range(len(pre_counts)):
        stats.append([f"{bin_edges[i]}-{bin_edges[i+1]}", pre_counts[i]])

    stats.append(["Mature Length", "Count"])
    for x, y in zip(mature_x, mature_y):
        stats.append([x, y])

    pd.DataFrame(stats).to_csv("miRNA_statistics.csv", index=False, header=False, encoding="utf-8-sig")

def main():
    args = parse_args()
    pre_seq, mature_seq, star_seq = load_data(args.input, CUSTOM_STYLE)
    bin_edges, pre_counts = get_pre_hist_data(pre_seq, CUSTOM_STYLE)
    cum_x, cum_y = get_cumulative_curve(pre_seq)
    mature_x, mature_y = get_mature_hist_data(mature_seq)

    plot_pre_hist_with_curve(bin_edges, pre_counts, cum_x, cum_y, CUSTOM_STYLE)
    plot_mature_hist(mature_x, mature_y, CUSTOM_STYLE)
    save_statistics(pre_seq, mature_seq, star_seq, bin_edges, pre_counts, mature_x, mature_y)
    print("down!")

if __name__ == "__main__":
    main()