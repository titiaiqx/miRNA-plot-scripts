# -*- coding: utf-8 -*-
"""
miRNA count 表达量 PCA 分析
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

FIG_SIZE = (9, 7)

GROUP_NAMES = ['FL', 'FS', 'ML', 'MS']

GROUP_COLORS = ['#E64B35', '#E64B35', '#4DBBD5', '#4DBBD5']

GROUP_MARKERS = ['o', '^', 'o', '^']

POINT_SIZE = 400

POINT_EDGE = True
POINT_EDGE_WIDTH = 1.2
POINT_EDGE_COLOR = 'black'

AXIS_LINE_WIDTH = 1.5
AXIS_LABEL_SIZE = 18
AXIS_TICK_SIZE = 16
AXIS_TITLE_SIZE = 16

SHOW_GRID = True
GRID_LINE_WIDTH = 0.6
GRID_STYLE = '--'

LEGEND_FONT_SIZE = 12
LEGEND_TITLE_SIZE = 13
LEGEND_MARKER_SIZE = 10

OUTPUT_PNG = "miRNA_PCA.png"
OUTPUT_SVG = "miRNA_PCA.svg"  
OUTPUT_CSV = "miRNA_PCA_results.csv"

parser = argparse.ArgumentParser(description='miRNA PCA Analysis')
parser.add_argument('-i', '--input', required=True, help='输入TPM CSV文件')
args = parser.parse_args()

df = pd.read_csv(args.input, index_col=0)
df_log = np.log2(df + 1)  

X = df_log.T.values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# PCA
pca = PCA(n_components=2)
pca_res = pca.fit_transform(X_scaled)
pc1_ratio = pca.explained_variance_ratio_[0] * 100
pc2_ratio = pca.explained_variance_ratio_[1] * 100

pca_df = pd.DataFrame({
    'Sample': df.columns,
    'PC1': pca_res[:, 0],
    'PC2': pca_res[:, 1]
})

def get_group(name):
    if name.startswith('FL'): return 'FL'
    elif name.startswith('FS'): return 'FS'
    elif name.startswith('ML'): return 'ML'
    elif name.startswith('MS'): return 'MS'
    else: return 'Other'
pca_df['Group'] = pca_df['Sample'].apply(get_group)

plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.figure(figsize=FIG_SIZE, dpi=300)

for i, group in enumerate(GROUP_NAMES):
    sub = pca_df[pca_df['Group'] == group]
    plt.scatter(
        sub['PC1'], sub['PC2'],
        s=POINT_SIZE,
        c=GROUP_COLORS[i],
        marker=GROUP_MARKERS[i],
        label=group,
        edgecolors=POINT_EDGE_COLOR if POINT_EDGE else None,
        linewidth=POINT_EDGE_WIDTH if POINT_EDGE else 0
    )

# 坐标轴设置
plt.xlabel(f'PC1 ({pc1_ratio:.2f}%)', fontsize=AXIS_LABEL_SIZE, weight='bold')
plt.ylabel(f'PC2 ({pc2_ratio:.2f}%)', fontsize=AXIS_LABEL_SIZE, weight='bold')
plt.title('miRNA Expression PCA Analysis', fontsize=AXIS_TITLE_SIZE, weight='bold', pad=15)

# 坐标轴粗细
ax = plt.gca()
ax.spines['top'].set_linewidth(AXIS_LINE_WIDTH)
ax.spines['bottom'].set_linewidth(AXIS_LINE_WIDTH)
ax.spines['left'].set_linewidth(AXIS_LINE_WIDTH)
ax.spines['right'].set_linewidth(AXIS_LINE_WIDTH)
plt.tick_params(axis='both', labelsize=AXIS_TICK_SIZE, width=AXIS_LINE_WIDTH)

if SHOW_GRID:
    plt.grid(alpha=0.3, linestyle=GRID_STYLE, linewidth=GRID_LINE_WIDTH)

plt.legend(
    title='Groups',
    loc='upper left',
    bbox_to_anchor=(1.02, 1),
    fontsize=LEGEND_FONT_SIZE,
    title_fontsize=LEGEND_TITLE_SIZE,
    frameon=True,
    shadow=False,
    labelspacing=1.2,
    handletextpad=0.8,
    borderpad=0.8
)
plt.setp(ax.get_legend().get_texts(), weight='bold')

plt.tight_layout()
plt.savefig(OUTPUT_SVG, format='svg', bbox_inches='tight')
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')
plt.close()

pca_df.to_csv(OUTPUT_CSV, index=False, encoding='utf-8-sig')

print("="*60)
print("down!")
