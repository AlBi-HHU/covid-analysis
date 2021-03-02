import sys
sys.path.append("..") 
import itertools

import csv
from collections import Counter

from shared import rev_comp

def main(reference, reads, th_min, th_frmr, th_max, th_covr, output, delKmers):

    #Parse deletion kmers from the input file into a set
    invalidDeletionKmers = set()
    for kmer in open(delKmers,'r').read().splitlines():
        invalidDeletionKmers.add(kmer)

    kmer_counts = read_jellyfish(reads)

    # Remove kmer can't be trust based on abundance and forward reverse ratio
    for kmer in generate_untrust_kmer(kmer_counts, th_max, th_min, th_frmr,invalidDeletionKmers):
        if kmer in kmer_counts:
            del kmer_counts[kmer]
        if rev_comp(kmer) in kmer_counts:
            del kmer_counts[rev_comp(kmer)] 

    # Remove kmer can't be trust based on around coverage
    reads_kmer = set()
    for kmer in kmer_counts.keys():
        count = kmer_counts[kmer] + kmer_counts[rev_comp(kmer)]

        # If kmer is very abundante we trust them
        if count > th_max:
            reads_kmer.add(kmer)
            continue
        
        preds = [kmer_counts[pred_mer] + kmer_counts[rev_comp(pred_mer)] for pred_mer in around_kmer(kmer, -1) if pred_mer in kmer_counts or rev_comp(pred_mer) in kmer_counts]
        succs = [kmer_counts[succ_mer] + kmer_counts[rev_comp(succ_mer)] for succ_mer in around_kmer(kmer, 1) if succ_mer in kmer_counts or rev_comp(succ_mer) in kmer_counts]

        if preds and succs:
            coverage = min(max(preds), max(succs))
        else:
            # if we can't guess the coverage kmer is probably an error
            continue

        if count > (coverage * th_covr):
            reads_kmer.add(kmer)
            reads_kmer.add(rev_comp(kmer))
            
    ref_kmer = set(read_jellyfish(reference).keys())

    print("trust reads kmer {}".format(len(reads_kmer)))
    print("trust reference kmer {}".format(len(ref_kmer)))
    print("total trust kmer {}".format(len(reads_kmer.union(ref_kmer))))

    with open(output, "w") as fh:
        for kmer in reads_kmer.union(ref_kmer):
            print(">1\n{}".format(kmer), file=fh)


def read_jellyfish(path):
    data = Counter()

    with open(str(path)) as fh:
        reader = csv.reader(fh, delimiter=' ')
        for row in reader:
            data[row[0]] = int(row[1])

    return data

def isHomopolymer(kmer,length=4):
    currentRepeatLength = 0
    lastBase = 'N'
    for base in kmer:
        if base == lastBase:
            currentRepeatLength += 1
            if currentRepeatLength == length:
                return True
        else:
            currentRepeatLength = 0
        lastBase = base
    return False

def generate_untrust_kmer(reads_kmer, th_max, th_min, th_frmr, invalidDeletions):
    """ generate kmer can't be trust based on abundance, homopolymer content and forward reverse ratio """

    for kmer in list(reads_kmer.keys()):
        if isHomopolymer(kmer):
            yield kmer
        else:
            forward = reads_kmer[kmer]
            reverse = reads_kmer[rev_comp(kmer)]

            total = forward + reverse

            if total > 0:
                ratio = min(forward / total, reverse / total)
            else:
                ratio = 0

            if total < th_min:
                yield kmer
            elif ratio < th_frmr:
                yield kmer
            else:
                for delKmer in invalidDeletions:
                    if delKmer in kmer:
                        yield kmer



def around_kmer(kmer, pos):
    if pos < 0:
        suffix = kmer[:pos]
        for prefix in generate_all_seq(abs(pos)):
            new_kmer = prefix + suffix
            if new_kmer != kmer:
                yield new_kmer
            
    elif pos > 0:
        prefix = kmer[pos:]
        for suffix in generate_all_seq(pos):
            new_kmer = prefix + suffix
            if new_kmer != kmer:
                yield new_kmer
            
    else:
        yield kmer


def generate_all_seq(length):
    for p in itertools.product(['A', 'C', 'T', 'G'], repeat=length):
        yield ''.join(p)


if "snakemake" in locals():
    th_min = int(config["kmerFilterMin"]),
    th_frmr = float(config["kmerFilterFRMinRatio"]),
    th_max = int(config["kmerFilterTrustAbundance"]),
    th_covr = float(config["kmerFilterGuessCoverageRatio"]),
    main(snakemake.input["reference"],
         snakemake.input["reads"],
         th_min,
         th_frmr,
         th_max,
         th_covr,
         snakemake.output["kmerset"],
         snakemake.input['delKmers'])
else:
    import sys

    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), sys.argv[7],sys.argv[8])
