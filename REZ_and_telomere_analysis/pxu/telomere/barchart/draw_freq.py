#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description='Draw rNMP freqs in telomere')
    parser.add_argument('freq', type=argparse.FileType('r'), help='rNMP incorporation frequency in telomere region')
    parser.add_argument('-o', help='Output figure name')
    args = parser.parse_args()

    # process data
    df = pd.read_csv(args.freq, sep='\t')
    g_columns = ['1_T', '2_T', '3_A','4_G', '5_G', '6_G']
    c_columns = ['-1_C', '-2_C', '-3_C','-4_T', '-5_A', '-6_A']
    data = []
    for _, r in df.iterrows():
        for c in g_columns:
            if r.Total_G == 0:
                continue
            name = c.split('_')[0]
            data.append([r['Sample'], r['Celltype'], name, r[c]/r['Total_G'], 'G'])
        for c in c_columns:
            if r.Total_C == 0:
                continue
            name = str(7+int(c.split('_')[0]))
            data.append([r['Sample'], r['Celltype'], name, r[c]/r['Total_C'], 'C'])
    data = pd.DataFrame(data, columns=['Sample', 'Celltype', 'Location', 'Frequency', 'Strand'])

    # Draw heatmap
    fig, ax = plt.subplots(nrows=2, figsize=(9, 6), dpi=300)
    sns.set(style='ticks', font_scale=8)
    plt.subplots_adjust(hspace=0.25)
    # build palette
    pal = ['#E31A1C', '#1F78B4', '#FFFFB9', '#33A02C']
    palette_g = {}
    for b in g_columns:
        i, base = b.split('_')
        color = pal['ACGT'.index(base)]
        palette_g[i] = color
    palette_c = {}
    for b in c_columns:
        i, base = b.split('_')
        color = pal['ACGT'.index(base)]
        palette_c[str(7+int(i))] = color
    # Barplot
    bars = sns.barplot(
        x='Celltype', 
        y='Frequency', 
        hue='Location', 
        hue_order=['1','2','3','4','5','6'],
        ax=ax[0],
        data=data[data.Strand=='G'],
        edgecolor='k',
        palette=palette_g,
        errwidth=1,
        capsize=0.1
        )
    # add sequence labels
    for i, label in enumerate('TTAGGG'):
        for p in bars.containers[i].patches:
            x = p.get_x()
            ax[0].text(x, -0.1, label, fontsize=20)
    bars = sns.barplot(
        x='Celltype', 
        y='Frequency', 
        hue='Location', 
        hue_order=['1','2','3','4','5','6'],
        ax=ax[1],
        data=data[data.Strand=='C'],
        edgecolor='k',
        palette=palette_c,
        errwidth=1,
        capsize=0.1
    )
    # add sequence labels
    for i, label in enumerate('AATCCC'):
        for p in bars.containers[i].patches:
            x = p.get_x()
            ax[1].text(x, 0, label, fontsize=20)
    # add data points
    sns.swarmplot(
        x='Celltype', 
        y='Frequency', 
        hue='Location', 
        hue_order=['1','2','3','4','5','6'],
        ax=ax[0],
        data=data[data.Strand=='G'],
        palette='dark:black',
        size=3,
        dodge=True
        )
    sns.swarmplot(
        x='Celltype', 
        y='Frequency', 
        hue='Location', 
        hue_order=['1','2','3','4','5','6'],
        ax=ax[1],
        data=data[data.Strand=='C'],
        palette='dark:black',
        size=3,
        dodge=True
    )
    # y scale
    ax[0].set_ylim([0,0.8])
    ax[1].set_ylim([0,0.8])
    # Y-axis of the bottom
    ax[1].invert_yaxis()
    # despine
    sns.despine(ax=ax[0])
    sns.despine(ax=ax[1], top=False, bottom=True)
    # legends and ticks
    ax[0].set_xlabel('')
    ax[1].set_xlabel('')
    ax[0].set_ylabel('')
    ax[1].set_ylabel('')
    ax[0].set_xticklabels('')
    ax[1].set_xticklabels([])
    # remove legend
    ax[0].get_legend().remove()
    ax[1].get_legend().remove()
    # ytick lable size
    ax[0].tick_params(axis='y', labelsize=16)
    ax[1].tick_params(axis='y', labelsize=16)
    fig.savefig(args.o)


    print('Done!')
if __name__ == '__main__':
    main()