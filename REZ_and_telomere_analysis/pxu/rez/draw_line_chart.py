#!/usr/bin/env python3

from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FixedLocator

# get chromosome sizes
def get_size(ef):
    data = ef.groupby('Chromosome').End.max().copy()
    data = data.reset_index()[['Chromosome', 'End']].set_index('Chromosome')['End']
    return data.to_dict()


# get rNMP enriched zones
def get_rez(ef, nbins, threshold, sample_threshold):
    # add celltype
    if 'Celltype' not in ef.columns:
        ef['Celltype'] = ef['Sample'].apply(sep_sample)
    # get bin size
    ef['length'] = ef['End'] - ef['Start']
    bin = ef.length.median()
    # calculated EF for each bin
    total_count = ef.groupby('Sample').Count.sum().to_dict()
    total_size = ef.groupby('Sample').length.sum().to_dict()
    # merge
    ef['Zone_start'] = ef.Start // (bin*nbins) * (bin*nbins)
    rezs = ef.groupby(
        by=['Sample', 'Celltype', 'Chromosome', 'Strand', 'Zone_start']
        ).agg({
            'End':'max',
            'Count':'sum'
        }).reset_index().rename(columns={'Zone_start':'Start'})
    rezs['EF'] = rezs.apply(calc_ef, axis=1, total=total_count, genome=total_size)
    rezs = rezs[rezs.EF > threshold]
    # common rNMP enriched zones
    n_sample = len(rezs.Sample.unique())
    n_celltype = len(rezs.Celltype.unique())
    common = rezs.groupby(by=['Chromosome', 'Strand', 'Start', 'End']).agg(
        {
            'Sample':['nunique','unique'],
            'Celltype':'nunique',
            'EF':['mean','median']
        }
    ).reset_index()
    common.columns = ['Chromosome', 'Strand', 'Start', 'End', 'Sample_count', 'Sample', 'Celltype_count', 'EF_mean', 'EF_median']
    common = common[(common.Celltype_count == n_celltype) & (common.Sample_count > n_sample * sample_threshold)].copy()
    common = common.sort_values(by=['Strand', 'Start'])
    names = []
    curr = defaultdict(lambda: {'+':1, '-':1})
    for _, x in common.iterrows():
        st = x.Strand
        chrom = x.Chromosome
        if st == '+':
            names.append(f'REZ{chrom}-{curr[chrom][st]}F')
        else:
            names.append(f'REZ{chrom}-{curr[chrom][st]}R')
        curr[chrom][st] += 1
    common['Name'] = names
    common = common[['Name', 'Chromosome', 'Strand', 'Start', 'End', 'Sample_count', 'EF_mean', 'EF_median', 'Sample']]
    common['Sample'] = common['Sample'].apply(concat_sample)
    common = common.sort_values(by=['Strand', 'Name'])
    return rezs, common


# concat sample
def concat_sample(x):
    return ','.join(x)

# add sample
def sep_sample(x):
    return x.split('_')[0]

# calculate ef
def calc_ef(x, total, genome):
    return (x.Count/(x.End - x.Start))/(total[x.Sample]/genome[x.Sample])


# draw figures
def draw(ef, common_rezs, size, output, tick_interval=20000000):
    # add celltype
    if 'Celltype' not in ef.columns:
        ef['Celltype'] = ef['Sample'].apply(sep_sample)
    # get bin size
    ef['length'] = ef['End'] - ef['Start']
    bin = ef.length.median()
    # calculated EF for each bin
    total_count = ef.groupby('Sample').Count.sum().to_dict()
    total_size = ef.groupby('Sample').length.sum().to_dict()
    ef['EF'] = ef.apply(calc_ef, axis=1, total=total_count, genome=total_size)
    # Draw figures
    nchroms = len(size)
    sns.set(style='ticks', font_scale=1.8)
    fig, ax = plt.subplots(nrows=nchroms, figsize=(10, 15), dpi=600)
    plt.subplots_adjust(left=0.15, right=0.98, top=0.98, bottom=0.03, hspace=0.3)
    # define palette
    palette = {
        'CD4T': "#ff7f0e",
        'hESC': "#2ca02c",
        # 'WB-GTP Control':"#ff7f0e",
        # 'WB-GTP PTSD': "#d62728",
        # 'HCT116': "#2ca02c",
        'HEK293T': "#d62728",
        'RNH2A-KO T3-8': "#e377c2",
        'RNH2A-KO T3-17': "#9467bd"
    }
    # add other types
    color_cands = sns.color_palette('Set2')
    i = 0
    for k in ef.Celltype.unique():
        if k not in palette:
            palette[k] = color_cands[i]
            i += 1
    efs = ef.groupby(['Chromosome', 'Strand'])
    commons = common_rezs.groupby(['Chromosome', 'Strand'])
    xlim = max(size.values())+1000000
    for i, chrom in enumerate(sorted(size.keys(), key=lambda x:int(x.replace('chr', '')))):
        print(f'Start draw {chrom}')
        ef_fwd = efs.get_group((chrom, '+'))
        ef_rev = efs.get_group((chrom, '-')).copy()
        ef_rev['EF'] = -ef_rev['EF']
        sns.lineplot(x='Start', y='EF', hue='Celltype', data=ef_fwd, 
            palette=palette, ax=ax[i], legend=None, linewidth=0.3)
        sns.lineplot(x='Start', y='EF', hue='Celltype', data=ef_rev,
            palette=palette, ax=ax[i], legend=None, linewidth=0.3)
        # highlight common rezs
        if (chrom, '+') in commons.groups:
            common_rezs_fwd = commons.get_group((chrom, '+'))
            for _, x in common_rezs_fwd.iterrows():
                ax[i].axvspan(x.Start, x.End, color='yellow', alpha=0.7, ymin=0.5)
        if (chrom, '-') in commons.groups:
            common_rezs_rev = commons.get_group((chrom, '-'))
            for _, x in common_rezs_rev.iterrows():
                ax[i].axvspan(x.Start, x.End, color='yellow', alpha=0.7, ymax=0.5)
        # line for chrom sizes
        ax[i].plot([0, size[chrom]], [0, 0], 'k-', linewidth=1.5)
        ax[i].plot([size[chrom], size[chrom]], [-5, 5], 'k-', linewidth=1.5)
        # format
        sns.despine(ax=ax[i], bottom=True)
        ax[i].set_ylim((-4.5, 4.5))
        ax[i].yaxis.set_major_locator(FixedLocator([-4, 0, 4]))
        ax[i].set_yticklabels([4, 0, 4])
        ax[i].set_xlim((0, xlim))
        ax[i].xaxis.set_visible(False)
        ax[i].set_xlabel('')
        ax[i].set_ylabel(chrom, rotation='horizontal')
        ax[i].yaxis.set_label_coords(-0.12, -0.0)
    # overall formatting
    ax[-1].xaxis.set_visible(True)
    ax[-1].xaxis.set_major_locator(FixedLocator(list(range(0, xlim, tick_interval))))
    ax[-1].set_xticklabels([])
    fig.savefig(f'{output}_rezs.png')


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Draw REZ linechart from binned counts')
    parser.add_argument('bin_counts', type=argparse.FileType('r'), help='Bin counts TSV file')
    parser.add_argument('-o', default='REZ', help='Output basename, default=REZ')
    parser.add_argument('-n', type=int, default=10, help='# of bins for REZ size, default= 10 bins')
    parser.add_argument('--ef_threshold', type=float, default=1.5, help='Enrichment factor threshold for enriched regions, default=1.5')
    parser.add_argument('--sample_threshold', type=float, default=0.8, help='The minimum library ratio threshold of common enriched regions, default=0.8')
    parser.add_argument('--selected', default=None, nargs='+', help='Selected paricular cell type(s)')
    args = parser.parse_args()

    # load data
    ef = pd.read_csv(args.bin_counts, sep='\t')

    # get chromosome sizes
    size = get_size(ef)

    # enriched zones
    rezs, common = get_rez(ef, args.n, args.ef_threshold, args.sample_threshold)
    rezs.to_csv(f'{args.o}_rezs.tsv', sep='\t', index=False)
    common.to_csv(f'{args.o}_common_rezs.tsv', sep='\t', index=False)

    # draw figures
    draw(ef, common, size, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
