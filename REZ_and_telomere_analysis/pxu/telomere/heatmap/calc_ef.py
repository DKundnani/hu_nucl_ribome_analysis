#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Normalization of rNMP preference analysis')
    parser.add_argument('raw', type=argparse.FileType('r'), help='rNMP incorporation file for desired chromosome.')
    parser.add_argument('bg', type=argparse.FileType('r'), help='Background frequency')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file.')
    parser.add_argument('--name', default='saccer', help='Name of chromosome in background frequency used for normalization, default = saccer')
    args = parser.parse_args()

    name = args.bg.readline().rstrip('\n').split()[1:]

    # bg[species][di/tri] = freq
    bg = defaultdict(dict)
    for l in args.bg:
        ws = l.rstrip('\n').split()
        for i in range(1,len(ws)):
            bg[ws[0]][name[i-1]] = float(ws[i])

    # header
    l = args.raw.readline()
    args.o.write(l)
    di = l.rstrip('\n').split('\t')

    # normalize
    for l in args.raw:
        ws = l.split('\t')
        if len(ws) != len(di):
            continue
        # normalize freq
        freq = []
        sp = args.name
        assert sp in bg, f'[ERROR] Cannot find chrom {sp} in background file'
        # calc ef
        for i in range(1, len(ws)):
            if bg[sp][di[i]] == 0:
                freq.append(0)
            else:
                freq.append(float(ws[i])/bg[sp][di[i]])
        args.o.write(ws[0] + '\t' + '\t'.join([str(freq[i]/sum(freq)) for i in range(len(freq))]) + '\n')
    print('Done!')


if __name__ == '__main__':
    main()
