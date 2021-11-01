#!/usr/bin/env python 
from collections import Counter
import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-i',
                '--infile',
                metavar='File',
                help='Input bamfile',
                type=str,
                required=True)

parser.add_argument('-o',
                '--outfile',
                metavar='File',
                help='Output bamfile',
                type=str,
                required=True)

parser.add_argument('-f', 
                    '--frac',
                    help='Fraction threshold of clipped nucleotides above which to remove read',
                    default=0.5,
                    type=float,
                   required = False)

def main():
    with pysam.AlignmentFile(args.infile, "rb") as infile, pysam.AlignmentFile(args.outfile, "wb", header=infile.header) as outfile:
        for read in infile.fetch():
            c = Counter()
            for i,j in read.cigartuples:
                c[i] += j
            #S   BAM_CSOFT_CLIP  4
            #H   BAM_CHARD_CLIP  5
            if (c[4]+c[5])/sum(c.values()) < args.frac:
                _ =outfile.write(read)

args = parser.parse_args()
if __name__ == "__main__":
   main()