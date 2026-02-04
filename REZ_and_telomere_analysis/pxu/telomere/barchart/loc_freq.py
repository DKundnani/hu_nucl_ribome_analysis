#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

# read lib info
def read_info(fr):
    data = defaultdict(str)
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        data[ws[1]] = ws[0]
    return data


def rc(s):
    bases = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    res = ''
    for c in s:
        res = bases[c] + res
    return res


# Count rNMPs at each location of telomere
def count_rNMPs(data, fq, umi_length):
    lib, st = fq.name.split('/')[-1].split('.')[0].split('_')
    if st == 'c':
        seq = 'CCCTAA'
    else:
        seq = 'TTAGGG'
    if lib not in data:
        data[lib] = {'TTAGGG':[set() for _ in range(6)], 'CCCTAA':[set() for _ in range(6)]}
    # check each sequence
    cnt = 0
    for l in fq:
        cnt += 1
        if cnt % 4 == 1 and umi_length == 0:
            umi = l.rstrip('\n').split('_')[-1]
        if cnt % 4 != 2:
            continue
        # deduplicate with UMI
        if umi_length != 0:
            umi = l[:umi_length]
            l = l.rstrip('\n')[umi_length:]
            for i in range(6):
                if seq[i:] == l[:6-i]:
                    data[lib][rc(seq)][5-i].add(umi)
                    break
        else:
            for i in range(6):
                if seq[i:] == l[:6-i]:
                    data[lib][rc(seq)][5-i].add(umi)
                    break


def main():
    parser = argparse.ArgumentParser(description='Calculate the count of rNMPs at each location of the telomere')
    parser.add_argument('fq', nargs='+', type=argparse.FileType('r'), help='Input rNMP fastq files')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Order file containing library information')
    parser.add_argument('-o', default=sys.stdout, type=argparse.FileType('w'), help='Output to file')
    parser.add_argument('--umi_length', default=11, type=int, help='UMI+Barcode length, if set to 0, fetch umi from the read rame. default=11')
    args = parser.parse_args()

    # count freq
    data = {}
    for fq in args.fq:
        count_rNMPs(data, fq, args.umi_length)

    # libinfo
    libinfo = read_info(args.libinfo)

    # output to file
    args.o.write('Sample\tCelltype\tTotal\tTotal_G\t1_T\t2_T\t3_A\t4_G\t5_G\t6_G\tTotal_C\t-1_C\t-2_C\t-3_C\t-4_T\t-5_A\t-6_A\n')
    for lib, v in data.items():
        args.o.write(f'{lib}\t{libinfo[lib]}\t')
        total_g = sum([len(x) for x in v['TTAGGG']])
        total_c = sum([len(x) for x in v['CCCTAA']])
        total = total_g + total_c
        args.o.write(f'{total}\t{total_g}')
        for uniq in v['TTAGGG']:
            args.o.write(f'\t{len(uniq)}')
        args.o.write(f'\t{total_c}')
        for uniq in v['CCCTAA']:
            args.o.write(f'\t{len(uniq)}')
        args.o.write('\n')


    print('Done!')

if __name__ == '__main__':
    main()