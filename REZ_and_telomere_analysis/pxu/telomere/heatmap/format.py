#!/usr/bin/env python3

import argparse
from collections import defaultdict


def main():
    # argparse
    parser = argparse.ArgumentParser(description='Format frequencies for telomere')
    parser.add_argument('raw', type=argparse.FileType('r'), help='rNMP frequency in telomere.')
    parser.add_argument('-o', default='telomere', help='Output file base name')
    args = parser.parse_args()

    header = args.raw.readline().rstrip('\n').split()

    with open(args.o + '_mono_both.raw', 'w') as mono, \
        open(args.o + '_dinuc_d1_nr_both.raw', 'w') as dinuc, \
        open(args.o + '_mono_cstrand.raw', 'w') as mono_c, \
        open(args.o + '_dinuc_d1_nr_cstrand.raw', 'w') as dinuc_c, \
        open(args.o + '_mono_gstrand.raw', 'w') as mono_g, \
        open(args.o + '_dinuc_d1_nr_gstrand.raw', 'w') as dinuc_g:
        ms = list('ACGT')
        ds = ['AA', 'CA', 'GA', 'TA', 'AC', 'CC', 'GC', 'TC',
                'AG', 'CG', 'GG', 'TG', 'AT', 'CT', 'GT', 'TT']
        # header
        mono.write('Sample\t' + '\t'.join(ms) + '\n')
        dinuc.write('Sample\t' + '\t'.join(ds) + '\n')
        mono_c.write('Sample\t' + '\t'.join(ms) + '\n')
        dinuc_c.write('Sample\t' + '\t'.join(ds) + '\n')
        mono_g.write('Sample\t' + '\t'.join(ms) + '\n')
        dinuc_g.write('Sample\t' + '\t'.join(ds) + '\n')
        for l in args.raw:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 10:
                continue
            name = ws[0]
            data = defaultdict(int)
            data_c = defaultdict(int)
            data_g = defaultdict(int)
            # iterate each column
            gs = 'GTTAGGG'
            prefix_g = 3
            for i, x in enumerate(gs):
                if not i:
                    continue
                data[gs[i]] += int(ws[i+prefix_g])
                data[gs[i-1:i+1]] += int(ws[i+prefix_g])
                data_g[gs[i]] += int(ws[i+prefix_g])
                data_g[gs[i-1:i+1]] += int(ws[i+prefix_g])
            cs = 'ACCCTAA'
            prefix_c = 10
            for i, x in enumerate(cs):
                if not i:
                    continue
                data[cs[i]] += int(ws[i+prefix_c])
                data[cs[i-1:i+1]] += int(ws[i+prefix_c])
                data_c[cs[i]] += int(ws[i+prefix_c])
                data_c[cs[i-1:i+1]] += int(ws[i+prefix_c])
            mono.write(f'{ws[0]}\t' + '\t'.join([str(data[x]) for x in ms]) + '\n')
            dinuc.write(f'{ws[0]}\t' + '\t'.join([str(data[x]) for x in ds]) + '\n')
            mono_c.write(f'{ws[0]}\t' + '\t'.join([str(data_c[x]) for x in ms]) + '\n')
            dinuc_c.write(f'{ws[0]}\t' + '\t'.join([str(data_c[x]) for x in ds]) + '\n')
            mono_g.write(f'{ws[0]}\t' + '\t'.join([str(data_g[x]) for x in ms]) + '\n')
            dinuc_g.write(f'{ws[0]}\t' + '\t'.join([str(data_g[x]) for x in ds]) + '\n')
            

    print('Done!')


if __name__ == '__main__':
    main()
