# -*- coding: utf-8 -*-
"""
miRNA count 样本相关性热图 + 层次聚类
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.collections as mcoll
from matplotlib.colors import LinearSegmentedColormap
import warnings
import os
warnings.filterwarnings('ignore')

COLOR_ANCHORS = [
    (0.48, "#055282"),
    (0.55, "#4C9AC9"),
    (0.75, "#ffffff"),
    (0.94, "#E2745E"),
    (1.00, "#841E0A")   
]

V_MIN = 0.48
V_MAX = 1.0

FIG_SIZE = (11.3, 10)

SHOW_VALUE = True
VALUE_SIZE = 18
VALUE_WEIGHT = "normal"

LABEL_SIZE = 19
LABEL_WEIGHT = "normal"
ROTATE_X = 0

TICK_LINE_WIDTH = 1.5

BORDER_WIDTH = 0
BORDER_COLOR = "white"

SHOW_CLUSTER = True
CLUSTER_LINE_WIDTH = 2

SHOW_COLOR_BAR = True
COLOR_BAR_TITLE = "Pearson Correlation"
COLOR_BAR_TITLE_FONT_WEIGHT = "normal"
COLOR_BAR_TITLE_FONT_SIZE = 15
COLOR_BAR_X = 1.07
COLOR_BAR_Y = 0.18
COLOR_BAR_WIDTH = 0.03
COLOR_BAR_HEIGHT = 0.5
COLOR_BAR_FONT_SIZE = 16
COLOR_BAR_FONT_WEIGHT = "normal"
COLOR_BAR_TICK_WIDTH = 1.5

PLOT_TITLE = "miRNA Expression Sample Correlation + Clustering"
TITLE_SIZE = 16
TITLE_WEIGHT = "bold"

SVG_FILE = "sample_cluster_heatmap.svg"
PNG_FILE = "sample_cluster_heatmap.png"

sorted_anchors = sorted(COLOR_ANCHORS, key=lambda x: x[0])
anchor_vals = [p[0] for p in sorted_anchors]
anchor_cols = [p[1] for p in sorted_anchors]
norm_pos = [(v - V_MIN) / (V_MAX - V_MIN) for v in anchor_vals]
color_steps = list(zip(norm_pos, anchor_cols))

custom_cmap = LinearSegmentedColormap.from_list(
    "custom_multi_anchor_cmap",
    color_steps,
    N=300
)

plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['svg.fonttype'] = 'none'

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True)
args = parser.parse_args()

df = pd.read_csv(args.input, index_col=0)
df_log = np.log2(df + 1)
corr = df_log.corr(method='pearson')

g = sns.clustermap(
    corr,
    cmap=custom_cmap,
    vmin=V_MIN,
    vmax=V_MAX,
    annot=SHOW_VALUE,
    fmt=".2f",
    annot_kws={"size": VALUE_SIZE, "weight": VALUE_WEIGHT},
    linewidths=BORDER_WIDTH,
    linecolor=BORDER_COLOR,
    cbar=SHOW_COLOR_BAR,
    row_cluster=SHOW_CLUSTER,
    col_cluster=SHOW_CLUSTER,
    dendrogram_ratio=(0.12, 0.12),
    figsize=FIG_SIZE
)

for col in g.ax_col_dendrogram.collections:
    if isinstance(col, mcoll.LineCollection):
        col.set_linewidth(CLUSTER_LINE_WIDTH)
for col in g.ax_row_dendrogram.collections:
    if isinstance(col, mcoll.LineCollection):
        col.set_linewidth(CLUSTER_LINE_WIDTH)

g.cax.set_position([COLOR_BAR_X, COLOR_BAR_Y, COLOR_BAR_WIDTH, COLOR_BAR_HEIGHT])
g.cax.tick_params(labelsize=COLOR_BAR_FONT_SIZE, width=COLOR_BAR_TICK_WIDTH, length=4)
for tick in g.cax.get_yticklabels():
    tick.set_fontweight(COLOR_BAR_FONT_WEIGHT)
g.cax.set_title(COLOR_BAR_TITLE, fontsize=COLOR_BAR_TITLE_FONT_SIZE, fontweight=COLOR_BAR_TITLE_FONT_WEIGHT, y=1.01)

g.ax_heatmap.tick_params(width=TICK_LINE_WIDTH, length=4)

g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    fontsize=LABEL_SIZE,
    fontweight=LABEL_WEIGHT,
    rotation=ROTATE_X,
    ha='center'
)

g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(),
    fontsize=LABEL_SIZE,
    fontweight=LABEL_WEIGHT,
    rotation=0,
    va='center'
)

g.fig.suptitle(
    PLOT_TITLE,
    fontsize=TITLE_SIZE,
    fontweight=TITLE_WEIGHT,
    x=0.5,
    y=1.05,
    ha='center'
)

plt.savefig(SVG_FILE, format="svg", bbox_inches="tight")
plt.savefig(PNG_FILE, dpi=300, bbox_inches="tight")
plt.close()

for val, col in sorted_anchors:
    print(f"数值 {val:.2f} → {col}")
print(f"down!")