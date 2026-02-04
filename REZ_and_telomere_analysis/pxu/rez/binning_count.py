#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

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
    for fr in bed:
        fs = fr.name.split('/')[-1].split('.')[0]
        if fs not in info:
            print(f'Library {fs} doesn\'t have related information, skipped')
            continue
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
            # only last bin have a different bin size
            if loc // bin == size[chrom] // bin:
                b = size[chrom] % bin
            else:
                b = bin
            data[fs][chrom][st][loc//bin] += 1
    df = []
    for fs, v in data.items():
        for chrom, v1 in v.items():
            for st, cnts in v1.items():
                for i, cnt in enumerate(cnts):
                    df.append([fs, info[fs], chrom, st, i*bin, min(i*bin+bin, size[chrom]), cnt])
    df = pd.DataFrame(df, columns=['Sample', 'Celltype', 'Chromosome', 'Strand', 'Start','End','Count'])
    return df


# read order
def read_libinfo(fr, c):
    data = {}
    for l in fr:
        ws = l.rstrip().split('\t')
        name = ws.pop(c)
        data[name] = '_'.join(ws)
    return data


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Count abundance of rNMPs in each bin')
    parser.add_argument('bed', type=argparse.FileType('r'), nargs='+', help='rNMP incorporation bed files')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Library information')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file.')
    parser.add_argument('-c', type=int, default=1, help='Col num for FS number, default=1')
    parser.add_argument('-b', type=int, default=10000, help='Bin size, default=10,000 nt')
    parser.add_argument('--selected', default=None, nargs='+', help='Selected paricular cell type(s)')
    parser.add_argument('--removed_chromosomes', default=['chrM', 'chrX', 'chrY'], type=list, help='Chromosome to be removed, default=chrX, chrY, chrM')
    args = parser.parse_args()
    args.c -= 1

    # mt size
    size = get_chrom_size(args.fai, args.removed_chromosomes)

    # load information
    info = read_libinfo(args.libinfo, args.c)

    # load data
    cnts = load_bed(args.bed, size, info, args.b, args.selected)
    cnts.to_csv(args.o, sep='\t', index=False)

    print('Done!')

if __name__ == '__main__':
    main()
