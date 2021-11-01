#!/usr/bin/env python 
from pysam import VariantFile
import argparse
import sys
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i',
                '--infile',
                metavar='File',
                help='Input bamfile; should be sorted and indexed',
                type=str,
                required=True)

parser.add_argument('-o',
                '--outfile',
                metavar='File',
                help='Output fasta file for consensus (default stdout)',
                default=sys.stdout,
                type=argparse.FileType('w'),
                required=False)

args = parser.parse_args()

if __name__ == '__main__':
    with args.outfile as outfile:
        print("run_accession","CHROM","POS","DP","REF","ALT","REF_fwd","REF_rev","ALT_fwd","ALT_rev","AF","MAJOR","INDEL", sep="\t", file=outfile)
        vcf_file = VariantFile(args.infile)

        for rec in vcf_file.fetch():
            if rec.alleles == None:
                continue
            ref_allele = rec.alleles[0]
            ref_fwd = rec.info['DP4'][0]
            ref_rev = rec.info['DP4'][1]

            major = rec.info['MAJOR']
            indel = rec.info['INDEL']

            af = rec.info['AF']

            for i, allele in enumerate(rec.alleles):
                #Skip the reference allele
                if i == 0:
                    continue

                #Allele count
                allele_fwd = rec.info['DP4'][i*2]
                allele_rev = rec.info['DP4'][(i*2)+1]
                
                print('\t'.join([os.path.basename(args.infile).split('.')[0], rec.chrom, str(rec.pos), str(rec.info["DP"]), ref_allele, allele, str(ref_fwd), str(ref_rev), str(allele_fwd), str(allele_rev), str((allele_fwd+allele_rev)/rec.info["DP"]), str(major), str(indel)]), file = outfile)
