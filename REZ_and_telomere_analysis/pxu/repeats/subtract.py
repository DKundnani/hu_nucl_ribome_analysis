#!/usr/bin/env python3

import argparse
import sys

def main():
    # argparse
    parser = argparse.ArgumentParser(description='Subtract two raw rNMP bed files')
    parser.add_argument('bed1', type=argparse.FileType('r'), help='Larger bed file')
    parser.add_argument('bed2', type=argparse.FileType('r'), help='Smaller bed file')
    parser.add_argument('-o', type=argparse.FileType('w'), default=sys.stdout, help='Output to file')
    parser.add_argument('-i', type=int, default=4, help='Column number of read id, (4)')
    args = parser.parse_args()
    args.i -= 1

    # read smaller bed file
    seen = set()
    for l in args.bed2:
        ws = l.rstrip('\n').split('\t') 
        if len(ws) < args.i+1:
            continue
        seen.add(ws[args.i])
    
    # deal with larger file
    for l in args.bed1:
        ws = l.rstrip('\n').split('\t') 
        if len(ws) < args.i+1:
            continue
        if not ws[args.i] in seen:
            args.o.write(l)
    
    print('Done!')
        

if __name__ == '__main__':
    main()
