#!/usr/bin/env python 
import argparse
import pysam

parser = argparse.ArgumentParser()

parser.add_argument('-i',
                '--infile',
                help='Input bamfile (default stdin)',
                type=str,
                required=True)

parser.add_argument('-o',
                '--outfile',
                help='Output bamfile (default stout)',
                type=str,
                required=True)

parser.add_argument('-b',
                '--bedfile',
                help='Input primer bedfile containing the columns chrom, start, end, names and strand',
                type=str,
                required=True)

parser.add_argument('-p', 
                    '--padding',
                    help='Upstream padding added to primer location when searching for overlapping reads',
                    default=0,
                    type=int,
                   required = False)

class Primer(object):
    """Primer location class"""
    def __init__(self, arg):
        bed_values = arg.strip().split('\t')
        try:
            self.chrom = bed_values[0]
            self.start = int(bed_values[1])
            self.end = int(bed_values[2])
            if self.start >= self.end:
                raise ValueError("Something went wrong: primer start should be smaller than primer end")
            self.len = self.end - self.start
            self.names = bed_values[3]
            if bed_values[4] == '+':
                self.strand = "left"
            elif bed_values[4] == '-':
                self.strand = "right"
            else:
                raise ValueError("Something went wrong: ["+bed_values[4]+"] is not a correct strand indicator, please use + or -")
        except IndexError:
            print("Something went wrong: In your bedfile ["+arg+"] does not contain chrom, start, end, names and strand which is needed for primer clipping!")

def calculate_overlap(read, primer, padding):
    if primer.strand == 'left':
        if read.get_overlap(primer.start, primer.end + 1) > 0:
            if (read.reference_start > (primer.start - padding)):
                n_clip = primer.end - (read.reference_start - 1)
                if n_clip < 0:
                    raise ValueError("Something went wrong: Trying to clip a negative number of nucleotides in read "+read.query_name)
                return(n_clip)
    elif primer.strand == 'right':
        if read.get_overlap(primer.start, primer.end) > 0:
            if (read.reference_end < (primer.end + padding)):
                n_clip = read.reference_end - primer.start
                if n_clip < 0:
                    raise ValueError("Something went wrong: Trying to clip a negative number of nucleotides in read "+read.query_name)
                return(n_clip)
    #Return 0 if no overlap or overlap is not at the ends of the read
    return(0)

def clip_read(read, n_clip, side):
    '''
    Pysam cigar codes:
    M   BAM_CMATCH  0
    I   BAM_CINS    1
    D   BAM_CDEL    2
    N   BAM_CREF_SKIP   3 (not handled)
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD    6 (not handled)
    =   BAM_CEQUAL  7 (not handled)
    X   BAM_CDIFF   8 (not handled)
    B   BAM_CBACK   9 (not handled)
    '''
    
    #If all aligned bases are clipped completely clip cigarstring
    if read.qlen <= n_clip:
        read.cigartuples = [(4,read.query_length)]
        return

    current_cigar = read.cigartuples
    #Reverse cigar if we have to clip the right side
    if side == "right":
        current_cigar.reverse()

    #Expand cigarstring
    cigar_expanded = ''.join([j*str(i) for i,j in current_cigar])

    #Clip cigar untill no more "clip" left:
    clip_leftover = n_clip
    n = 0
    cigar_replacement = ''
    while clip_leftover > 0:
        cig = int(cigar_expanded[n])
        if cig == 0:
            cigar_replacement += "4" #Replace matches with softclipped base
            clip_leftover -= 1
        elif cig == 1:
            cigar_replacement += "4" #Replace insertions with softclipped base
        elif cig == 2:
            n_clip += 1 #Increase n_clip to increase reference_start shift at the end
        elif cig == 4:
            cigar_replacement += "4" #Do not replace softclipped
        elif cig == 5:
            cigar_replacement += "5" #Do not replace hardclipped
        else:
            raise ValueError("Something went wrong: do not know how to clip " + str(cig) + " in cigarstring")
        n += 1

    cigar_expanded = cigar_replacement + cigar_expanded[n:]

    #Recreate tuples:
    clipped_cigar = []
    c = 0
    nprev=cigar_expanded[0]
    for n in cigar_expanded:
        if n==nprev:
            c+=1
        else:
            clipped_cigar.append((int(nprev),c))
            c=1
            nprev=n
    clipped_cigar.append((int(n),c))

    #Un-reverse cigar if clipped on right side
    if side == "right":
        clipped_cigar.reverse()

    #Replace cigar tuples of the read
    read.cigartuples = clipped_cigar

    #Shift alignment start if clipped to the left
    if side == 'left':
        read.reference_start = read.reference_start + n_clip

def main():
    primers = []
    with open(args.bedfile,"r") as bedfile:
        for n, line in enumerate(bedfile):
            #Skip the header if it's there
            if (n == 0) and ("chrom" in line.strip().split('\t')):
                continue
            primers.append(Primer(line))

    with pysam.AlignmentFile(args.infile, "rb") as infile, pysam.AlignmentFile(args.outfile, "wb", header=infile.header) as outfile:
        for read in infile.fetch():
            for primer in primers:
                if read.is_unmapped:
                    continue
                overlap = calculate_overlap(read, primer, padding=args.padding)
                if overlap > 0:
                    clip_read(read, overlap, primer.strand)
            _ =outfile.write(read)

args = parser.parse_args()
if __name__ == "__main__":
   main()