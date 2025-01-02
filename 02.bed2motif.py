import gzip
import bz2
import argparse
from collections import defaultdict
from itertools import zip_longest

parser = argparse.ArgumentParser(description='Analysis k-mer endmotif')
parser.add_argument('-b', '--bedfile', required=True, type=str, dest='bedfile', help='<str> Input bedfile.')
parser.add_argument('-o', '--out_prefix', required=True, type=str, dest='out_prefix', help='<str> Prefix of output files.')
parser.add_argument('-minL', '--minLen', required=False, type=int, default=20,   dest='minLen', help='<int> Minimum length of fragments. [1] bp.')
parser.add_argument('-maxL', '--maxLen', required=False, type=int, default=600, dest='maxLen', help='<int> Maximum length of fragments. [600] bp.')
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
args = (parser.parse_args())
print(args)


def read_file(file):
    if file.split('.')[-1] == 'bz2':
        fh = bz2.open(file, 'rt')
    elif file.split('.')[-1] == 'gz':
        fh = gzip.open(file, 'rt')
    else:
        fh = open(file, 'r')
    return fh


fraglen_cnt = defaultdict(int)  # statistics fragment length distribution
motif_cnt = defaultdict(int)  # statistics motif
banded = defaultdict(lambda: defaultdict(lambda: 0))  # ccca banded

fh = read_file(args.bedfile)
for line in fh:
    # chr1    10179871        10180049        178     0       +       TTTT    ACCG
    if len(line.strip().split())==8:
        chrom, frag_start, frag_end, frag_len, ot1, ot2, motif1, motif2 = line.strip().split()
        frag_start, frag_end, frag_len = map(int, [frag_start, frag_end, frag_len])
        assert frag_len == frag_end - frag_start
        if frag_len < args.minLen or frag_len > args.maxLen: continue
        fraglen_cnt[frag_len] += 1

        if len(set(motif1) - set('ATCG')) == 0:
            motif_cnt[motif1] += 1
            banded[frag_len][motif1] += 1

        if len(set(motif2) - set('ATCG')) == 0:
            motif_cnt[motif2] += 1
            banded[frag_len][motif2] += 1


# statistics fragment length
out2 = open('{}.size'.format(args.out_prefix), 'w')
out2.write('#Size\tCount\tPercent\n')
size_tot_cnt = sum(fraglen_cnt.values())
for frag_len in range(args.minLen, args.maxLen + 1):
#for frag_len in sorted(fraglen_cnt.keys()):
    out2.write('{}\t{}\t{:.4f}\n'.format(frag_len, fraglen_cnt[frag_len], fraglen_cnt[frag_len] / size_tot_cnt * 100))
out2.close()

# statistics motif frequency
out3 = open('{}.motif'.format(args.out_prefix), 'w')
out3.write("#Motif\tCount\tPercent\n")
motif_tot_cnt = sum(motif_cnt.values())
for motif in sorted(motif_cnt.keys()):
    out3.write('{}\t{}\t{:.4f}\n'.format(motif, motif_cnt[motif], motif_cnt[motif] / motif_tot_cnt * 100))
out3.close()

# statistics banded
out4 = open('{}.banded.ccca'.format(args.out_prefix), 'w')
out4.write('#Size\tMotif\tCount\tPercent\n')
for motif in sorted(motif_cnt.keys()):
    for frag_len in range(args.minLen, args.maxLen + 1):
        try:
            pct = banded[frag_len][motif] / fraglen_cnt[frag_len] * 50
        except ZeroDivisionError:
            pct = 0
            pass
        out4.write('{}\t{}\t{}\t{:.4f}\n'.format(frag_len, motif, banded[frag_len][motif], pct))  # one size havs two motifs
out4.close()

