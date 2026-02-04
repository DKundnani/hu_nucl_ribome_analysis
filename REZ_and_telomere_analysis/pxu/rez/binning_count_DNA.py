#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
from collections import defaultdict

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


# Calc sequencing count in each bin
def load_sam(sam, size, bin, max_insert_size):
    data = {}
    locs = {}
    for fr in sam:
        fs = fr.name.split('/')[-1].split('.')[0]
        # initiate locs
        if fs not in locs:
            locs[fs] = {
                k:{'+':defaultdict(int), '-':defaultdict(int)} \
                for k in size.keys()
            }
        # Iterations
        for l in fr:
            if l[0] == '#':
                continue
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 11:
                continue
            chrom = ws[2]
            if chrom not in size:
                continue
            # get positions
            s = int(ws[3]) - 1
            length = int(ws[8])
            if length <= 0 or length > max_insert_size:
                continue
            e = s + length
            # get strand
            flag = int(ws[1])
            st = 1
            # second in pair
            if (flag % 256)//128:
                st *= -1
            # reverse
            if (flag % 32)//16:
                st *= -1
            st = '+' if st > 0 else '-'
            # store start and end
            locs[fs][chrom][st][s] += 1
            locs[fs][chrom][st][e] -= 1
    # initiate output
    df = []
    for fs, v in locs.items():
        for chrom, v1 in v.items():
            for st, d in v1.items():
                d[size[chrom]] += 0
                coords = sorted(d.keys())
                curr = 0
                depth = 0
                cnt = 0
                for coord in coords:
                    # while not in the same bin, move curr to the start of next bin and output
                    while curr // bin != coord // bin:
                        new = min((curr // bin + 1) * bin, size[chrom])
                        cnt += (new - curr) * depth 
                        df.append([fs, chrom, st, curr//bin*bin, new, cnt])
                        # start new bin
                        curr = new
                        cnt = 0
                    # update cnt and depth
                    cnt += (coord - curr) * depth
                    depth += d[coord]
                    curr = coord
                # last bin
                df.append([fs, chrom, st, curr//bin*bin, size[chrom], cnt])
    df = pd.DataFrame(df, columns=['Sample', 'Chromosome', 'Strand', 'Start','End','Count'])
    return df


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Count coverage of DNA-seq sam file for each bin')
    parser.add_argument('sam', type=argparse.FileType('r'), nargs='+', help='Aligned_sam_files')
    parser.add_argument('fai', type=argparse.FileType('r'), help='Fasta index file')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file.')
    parser.add_argument('-b', type=int, default=10000, help='Bin size, default=10,000 nt')
    parser.add_argument('-m', type=int, default=1000, help='Max insert size allowed, default=1,000 nt')
    parser.add_argument('--removed_chromosomes', default=['chrM', 'chrX', 'chrY'], type=list, help='Chromosome to be removed, default=chrX, chrY, chrM')
    args = parser.parse_args()

    # mt size
    size = get_chrom_size(args.fai, args.removed_chromosomes)

    # load data
    cnts = load_sam(args.sam, size, args.b, args.m)
    cnts.to_csv(args.o, sep='\t', index=False)

    print('Done!')

if __name__ == '__main__':
    main()
