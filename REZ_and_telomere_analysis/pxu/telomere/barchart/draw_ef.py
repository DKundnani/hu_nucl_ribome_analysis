#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description='Draw telomere EF')
    parser.add_argument('ef', type=argparse.FileType('r'), help='Telomere EF CSV file')
    parser.add_argument('-o', help='Output figure basename')
    args = parser.parse_args()

    # process data
    df = pd.read_csv(args.ef)

    # Draw heatmap
    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
    sns.set(style='ticks', font_scale=8)
    # build palette
    pal = ['blue']
    # Barplot
    bars = sns.barplot(
        x='Celltype', 
        y='EF', 
        ax=ax,
        data=df[df.Strand!='Total'],
        edgecolor='k',
        palette=pal,
        errwidth=1,
        capsize=0.1
        )
    sns.swarmplot(
        x='Celltype', 
        y='EF', 
        ax=ax,
        data=df[df.Strand!='Total'],
        edgecolor='k',
        palette='dark:black',
        size=3
        )
    # plot a dashed line for EF=1
    xlims = ax.get_xlim()
    ax.plot(xlims, [1,1], 'k--')

    # despine
    sns.despine(ax=ax)
    # legends and ticks
    ax.set_xlabel('')
    ax.set_ylabel('')
    # ytick lable size
    ax.tick_params(axis='y', labelsize=32)
    ax.tick_params(axis='x', labelsize=0)
    fig.savefig(args.o + "_cgratio.png")


    # Draw heatmap
    fig, ax = plt.subplots(figsize=(9, 6), dpi=300)
    sns.set(style='ticks', font_scale=8)
    # build palette
    pal = ['darkorange']
    # Barplot
    bars = sns.barplot(
        x='Celltype', 
        y='EF', 
        ax=ax,
        data=df[df.Strand=='Total'],
        edgecolor='k',
        palette=pal,
        errwidth=1,
        capsize=0.1
        )
    sns.swarmplot(
        x='Celltype', 
        y='EF', 
        ax=ax,
        data=df[df.Strand=='Total'],
        edgecolor='k',
        palette='dark:black',
        size=3
        )
    
    # plot a dashed line for EF=1
    xlims = ax.get_xlim()
    ax.plot(xlims, [1,1], 'k--')
    # despine
    sns.despine(ax=ax)
    # legends and ticks
    ax.set_xlabel('')
    ax.set_ylabel('')
    # ytick lable size
    ax.tick_params(axis='y', labelsize=32)
    ax.tick_params(axis='x', labelsize=0)
    fig.savefig(args.o + "_both.png")

    print('Done!')
if __name__ == '__main__':
    main()