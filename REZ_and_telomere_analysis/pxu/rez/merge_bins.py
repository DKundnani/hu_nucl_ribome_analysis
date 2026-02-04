#!/usr/bin/env python3

import pandas as pd
import argparse
import sys


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Merge bins of counts')
    parser.add_argument('bin_counts', type=argparse.FileType('r'), help='Bin counts TSV file')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('-n', type=int, default=10, help='# of bins for merge, default= 10 bins')
    args = parser.parse_args()

    # load data
    cnts = pd.read_csv(args.bin_counts, sep='\t')

    # get bin size
    cnts['length'] = cnts['End'] - cnts['Start']
    bin = cnts.length.median()

    # merge
    cnts['Zone_start'] = (cnts.Start // (bin*args.n) * (bin*args.n)).astype(int)
    group = []
    for c in cnts.columns:
        if c not in ('End', 'Count', 'length', 'Start'):
            group.append(c)
    cnts = cnts.groupby(
        by=group
        ).agg({
            'End':'max',
            'Count':'sum'
        }).reset_index().rename(columns={'Zone_start':'Start'})

    # output
    cnts.to_csv(args.o, sep='\t', index=False)

    print('Done!')

if __name__ == '__main__':
    main()
