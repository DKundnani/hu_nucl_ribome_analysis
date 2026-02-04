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
def get_chrom_size(fai, names):
    # get mt length
    data = {}
    for l in fai:
        ws = l.split('\t')
        if ws[0] not in names:
            size = int(ws[1])
            data[ws[0]] = size
    return data


# load rNMP coordinates
def load_bed(bed, size, info, bin, selected):
    data = {}
    counts = {}
    total = sum(size.values())
    for fr in bed:
        fs = fr.name.split('/')[-1].split('.')[0]
        if fs not in info:
            print(f'Library {fs} doesn\'t have related information, skipped')
            continue
        counts[fs] = {'+':0, '-':0}
        geno = info[fs]
        if selected and geno not in selected:
            continue
        if fs not in data:
            data[fs] = {k:{} for k in size.keys()}
            for k in size.keys():
                data[fs][k] = {'+':[0]*(size[k]//bin+1), '-':[0]*(size[k]//bin+1)}
        for l in fr:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            chrom = ws[0]
            if chrom not in size:
                continue
            loc = int(ws[1])
            st = ws[5]
            counts[fs][st] += 1
            # only last bin have a different bin size
            if loc // bin == size[chrom] // bin:
                b = size[chrom] % bin
            else:
                b = bin
            data[fs][chrom][st][loc//bin] += 1/ (b/total)
    df = []
    for fs, v in data.items():
        for chrom, v1 in v.items():
            for st, efs in v1.items():
                for i, ef in enumerate(efs):
                    df.append([fs, info[fs], chrom, st, i*bin, min(i*bin+bin, size[chrom]), ef/counts[fs][st]])
    df = pd.DataFrame(df, columns=['Sample', 'Celltype', 'Chromosome', 'Strand', 'Start','End','EF'])
    return df


# read order
def read_libinfo(fr, c):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = '_'.join(ws)
    return data


# get rNMP enriched zones
def get_rez(ef, bin, nbins, threshold, sample_threshold):
    ef['Zone_start'] = ef.Start // (bin*nbins) * (bin*nbins)
    rezs = ef.groupby(
        by=['Sample', 'Celltype', 'Chromosome', 'Strand', 'Zone_start']
        ).agg({
            'End':'max',
            'EF':'mean'
        }).reset_index().rename(columns={'Zone_start':'Start'})
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


# draw figures
def draw(ef, common_rezs, size, output, tick_interval=20000000):
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
    parser = argparse.ArgumentParser(description='Analyze enriched zones for human nuclear DNA. MtDNA is excluded in the script')
    parser.add_argument('bed', type=argparse.FileType('r'), nargs='+', help='rNMP incorporation bed files')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Library information')
    parser.add_argument('-o', default='REZ', help='Output basename, default=REZ')
    parser.add_argument('-c', type=int, default=1, help='Col num for FS number, default=1')
    parser.add_argument('-b', type=int, default=1000000, help='Bin size, default=100,000 nt')
    parser.add_argument('-m', type=int, default=10, help='Minimum REZ size in number of bins, default= 10 bins')
    parser.add_argument('--ef_threshold', type=float, default=1.5, help='Enrichment factor threshold for enriched regions, default=1.5')
    parser.add_argument('--sample_threshold', type=float, default=0.8, help='The minimum library ratio threshold of common enriched regions, default=0.8')
    parser.add_argument('--selected', default=None, nargs='+', help='Selected paricular cell type(s)')
    parser.add_argument('--removed_chromosomes', default=['chrM', 'chrX', 'chrY'], type=list, help='Chromosome to be removed, default=chrX, chrY, chrM')
    args = parser.parse_args()
    args.c -= 1

    # mt size
    size = get_chrom_size(args.fai, args.removed_chromosomes)

    # load information
    info = read_libinfo(args.libinfo, args.c)

    # load data
    ef = load_bed(args.bed, size, info, args.b, args.selected)
    ef.to_csv(f'{args.o}_ef.tsv', sep='\t', index=False)

    # enriched zones
    rezs, common = get_rez(ef, args.b, args.m, args.ef_threshold, args.sample_threshold)
    rezs.to_csv(f'{args.o}_rezs.tsv', sep='\t', index=False)
    common.to_csv(f'{args.o}_common_rezs.tsv', sep='\t', index=False)
    
    # ef = pd.read_csv('REZ_500000_1_0.6_ef.tsv', sep='\t')
    # common = pd.read_csv('REZ_500000_1_0.6_common_rezs.tsv', sep='\t')

    # draw figures
    draw(ef, common, size, args.o)

    print('Done!')

if __name__ == '__main__':
    main()
