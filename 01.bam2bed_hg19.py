import os
import re
import sys
import gzip
import bz2
import time
import pysam
import argparse
from collections import defaultdict
from itertools import zip_longest
from pyfasta import Fasta
from pprint import pprint
import subprocess
#import multiprocessing
#from multiprocessing.dummy import Pool as ThreadPool
#from bx.intervals.intersection import Intersecter, Interval


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam', required=True, type=str, dest='bam', help='<str> Input bam file.')
parser.add_argument('-o', '--outprefix', required=True, type=str, dest='outprefix', help='<str> Prefix of output files.')
parser.add_argument('-m', '--max_mismatch', required=False, type=int, default=2, dest='max_mismatch', help='<str> Maximum mismatches allowed on a read. [2] bp.')
parser.add_argument('-minL', '--minLen', required=False, type=int, default=1,   dest='minLen', help='<int> Minimum length of fragments. [1] bp.')
parser.add_argument('-maxL', '--maxLen', required=False, type=int, default=600, dest='maxLen', help='<int> Maximum length of fragments. [600] bp.')
#parser.add_argument('-p', '--processes', required=False, type=int, default=2, dest='nproc', help='Number of processes to use in parallel.')
#parser.add_argument('-f', '--folder', required=False, type=str, default=os.getcwd(), dest='folder', help='Specific the folder to put in the result file.')
parser.add_argument('-q', '--quiet', required=False, action='store_true', dest='quiet', default=False, help="Don't print progress on stderr.")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
args = (parser.parse_args())


def disp(txt):
    if not args.quiet: sys.stderr.write('{}\t{}\n'.format(time.asctime(), txt))


def read_file(file):
    if file.split('.')[-1] == 'bz2':
        fh = bz2.open(file, 'rt')
    elif file.split('.')[-1] == 'gz':
        fh = gzip.open(file, 'rt')
    else:
        fh = open(file, 'r')
    return fh


def revcomp(seq):
    return seq.upper()[::-1].translate(str.maketrans('AGCT','TCGA'))


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):   # fetch() will only iterate over alignments in the SAM/BAM file
        if read.is_proper_pair:
            if (read.flag == 83) or (read.flag == 99) or (read.flag == 147)  or (read.flag == 163): 
                qname = read.query_name
                if qname not in read_dict:
                    if read.is_read1:
                        read_dict[qname][0] = read
                    else:
                        read_dict[qname][1] = read
                else:
                    if read.is_read1:
                        yield read, read_dict[qname][1]
                    else:
                        yield read_dict[qname][0], read
                    del read_dict[qname]


def cigar_align(cig_tup, in_list):
    aligned_list = []
    qIdx = 0  # index on query sequence
    for cig_num, count in cig_tup:  # [(0, 75), (4, 10)]
        if cig_num == 0:  # match
            aligned_list += in_list[qIdx: qIdx + count]
            qIdx += count
        elif cig_num == 1:  # insertion
            qIdx += count
        elif cig_num == 2:  # deletion
            aligned_list += list('N' * count)
        elif cig_num == 4:  # soft clipping
            qIdx += count

    return aligned_list


def get_frag(bamFile):
    '''
    Get fragment sequences from reads.
    '''
    samfile = pysam.Samfile(bamFile, "rb")
    for read1, read2 in read_pair_generator(samfile):

        # if read1.has_tag('XA') or read1.has_tag('SA'): continue
        # if read2.has_tag('XA') or read2.has_tag('SA'): continue
        if read1.mapping_quality < 10 or read2.mapping_quality < 10: continue
        cig_tup1, cig_tup2 = read1.cigar, read2.cigar
        if len(set([cig_num for cig_num, count in cig_tup1]) - set([0, 4])) != 0:  # No other cigar chars besides M, S
            # print(cig_tup1)
            continue
        if len(set([cig_num for cig_num, count in cig_tup2]) - set([0, 4])) != 0:  # No other cigar chars besides M, S
            # print(cig_tup2)
            continue
        # assert len(set([cig_num for cig_num, count in cig_tup1]) - set([0, 4])) == 0
        # I_cnt = sum([count for cig_num, count in cig_tup if cig_num == 1])
        # D_cnt = sum([count for cig_num, count in cig_tup if cig_num == 2])
        # S_cnt = sum([count for cig_num, count in cig_tup if cig_num == 4])
        # E_cnt = sum([count for cig_num, count in cig_tup if cig_num == 7])
        # X_cnt = sum([count for cig_num, count in cig_tup if cig_num == 8])
        # print(f'I: {I_cnt}\tD: {D_cnt}\tS: {S_cnt}\t=: {E_cnt}\tX: {X_cnt}')
        # fragLen = len(read1.query_sequence) - I_cnt - S_cnt + D_cnt
        S_cnt1 = sum([count for cig_num, count in cig_tup1 if cig_num == 4])
        S_cnt2 = sum([count for cig_num, count in cig_tup2 if cig_num == 4])

        chrom1, chrom2 = samfile.get_reference_name(read1.reference_id), samfile.get_reference_name(read2.reference_id)
        if chrom1 != chrom2: continue
        chrom = chrom1

        start1, start2 = read1.reference_start, read2.reference_start
        end1, end2 = start1 + read1.query_length - S_cnt1, start2 + read2.query_length - S_cnt2
        start = min(start1, start2)
        end = max(end1, end2)
        size = end - start
        if size < args.minLen or size > args.maxLen: continue

        strand1 = '-' if read1.is_reverse else '+'
        strand2 = '-' if read2.is_reverse else '+'
        if strand1 == '+' and strand2 == '-' and  start1 > start2: continue
        if strand1 == '-' and strand2 == '+' and  start1 < start2: continue

#         mis1, mis2 = read1.get_tag('NM'), read2.get_tag('NM')
#         if mis1 > args.max_mismatch or mis2 > args.max_mismatch: continue

        try:
            fragRef = ref_fa_dict[chrom][start: end].upper()
        except KeyError:
            print('{} not exists in ref. Skip this record.'.format(chrom))
            continue

        motif1 = fragRef[:4]
        # if len(set(motif1) - set('ATCG')) != 0: continue

        motif2 = revcomp(fragRef)[:4]
        # if len(set(motif2) - set('ATCG')) != 0: continue

        yield (chrom, start, end, size, motif1, motif2)


if __name__ == '__main__':
    ref_fa = 'hg19_HBV_EBV_lambda_pUC19.fa'  # default reference: hg19
    ref_fa_dict = Fasta(ref_fa)

    out = open(f'{args.outprefix}.bed', 'w')
    idx = 0
    for chrom, start, end, size, motif1, motif2 in get_frag(args.bam):
        idx += 1
        out.write(f'{chrom}\t{start}\t{end}\t{size}\t0\t+\t{motif1}\t{motif2}\n')

        if idx % 2000000 == 0:  # print progress
            disp(f'{idx} PE reads processed.')

    out.close()
    disp('Done.')

