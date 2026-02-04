#!/usr/bin/env python3

import argparse
import sys

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Extend bed file')
    parser.add_argument('bed', type=argparse.FileType('r'), help='rNMP bed file with sequence')
    parser.add_argument('-l', type=int, default=5, help='Length of both side to be extended (5)')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    args = parser.parse_args()
    
    for l in args.bed:
        ws = l.rstrip('\n').split('\t')
        if len(ws) < 6:
            continue
        ws[1] = max(int(ws[1]) - args.l, 0)
        ws[2] = int(ws[2]) + args.l
        args.o.write(f'{ws[0]}\t{ws[1]}\t{ws[2]}\t{ws[3]}\t{ws[4]}\t{ws[5]}\n')


if __name__ == '__main__':
    main()
