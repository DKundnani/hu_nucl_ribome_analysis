#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

# read info
def read_info(fr):
    data = {}
    for l in fr:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 2:
            continue
        data[ws[1]] = ws[0]
    return data


# check rNMP
def check(seq):
    seq = seq.upper()
    g = 'TTAGGG'
    c = 'CCCTAA'
    loc = seq.find(g)
    if loc != -1:
        i = 5 - loc
        return f'{i+1}_{g[i]}'
    loc = seq.find(c)
    if loc != -1:
        i = 5 - loc
        return f'-{i+1}_{c[i]}'
    return None
    

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Extend bed files')
    parser.add_argument('beds', type=argparse.FileType('r'), nargs='+', help='rNMP bed file(s)')
    parser.add_argument('libinfo', type=argparse.FileType('r'), help='Order file containing library information')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()

    # fetch lib info    
    info = read_info(args.libinfo)

    # read rNMPs
    rNMPs = {}
    for bed in args.beds:
        name = bed.name.split('/')[-1].split('.')[0]
        rNMPs[name] = defaultdict(int)
        for l in bed:
            ws = l.rstrip('\n').split('\t')
            if len(ws) < 6:
                continue
            seq = ws[4]
            loc = check(seq)
            if not loc:
                continue
            rNMPs[name][loc] += 1
        
    # header
    args.o.write('Sample\tCelltype\tTotal\tTotal_G\t1_T\t2_T\t3_A\t4_G\t5_G\t6_G\tTotal_C\t-1_C\t-2_C\t-3_C\t-4_T\t-5_A\t-6_A\n')
    for k in rNMPs.keys():
        total_c = sum([rNMPs[k][x] for x in rNMPs[k].keys() if x.startswith('-')])
        total = sum(rNMPs[k].values())
        args.o.write('\t'.join([
            k,
            info[k],
            str(total),
            str(total - total_c),
            str(rNMPs[k]['1_T']),
            str(rNMPs[k]['2_T']),
            str(rNMPs[k]['3_A']),
            str(rNMPs[k]['4_G']),
            str(rNMPs[k]['5_G']),
            str(rNMPs[k]['6_G']),
            str(total_c),
            str(rNMPs[k]['-1_C']),
            str(rNMPs[k]['-2_C']),
            str(rNMPs[k]['-3_C']),
            str(rNMPs[k]['-4_T']),
            str(rNMPs[k]['-5_A']),
            str(rNMPs[k]['-6_A'])
        ]) + '\n')

    print('Done!')


if __name__ == '__main__':
    main()
