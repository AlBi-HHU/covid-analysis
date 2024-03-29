import re
import csv
from collections import Counter
from math import log2

##### This contains shared functions or definitions that are used across multiple scripts

# The ambiguity letters for nucleotides
ambiguityLetters = {
    frozenset(("A")): "A",
    frozenset(("C")): "C",
    frozenset(("G")): "G",
    frozenset(("T")): "T",
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
    frozenset(("C", "G", "T")): "B",
    frozenset(("A", "G", "T")): "D",
    frozenset(("A", "C", "T")): "H",
    frozenset(("A", "C", "G")): "V",
    frozenset(("A", "C", "T", "G")): "N",
}

# An inversion of the map above to allow for quick reverse lookups
ambiguityLetters_inverted = {v: k for k, v in ambiguityLetters.items()}


def ambiguousBase(frozenset):
    if frozenset in ambiguityLetters:
        return ambiguityLetters[frozenset]
    else:
        return None


def isAmbiguous(base):
    if base in ["N", "A", "C", "G", "T"]:
        return False
    elif base in ambiguityLetters_inverted:
        return True
    return None


######## Graph Aligner Helper Functions
def get_node2seq(graph_path):
    node2seq = dict()

    with open(graph_path) as graph_fh:
        reader = csv.reader(graph_fh, delimiter="\t")
        for row in reader:
            if row[0] == "S" and row[1]:
                node2seq[row[1]] = row[2]

    return node2seq


def read_node2len(path):
    node2len = dict()

    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            node2len[row["node"]] = int(row["length"])

    return node2len


####### Various
def compute_entropy(seq):
    entropy = 0
    for k, v in Counter(seq).items():
        entropy += v / len(seq) * log2(v / len(seq))

    return entropy * -1


def strand_bias(record, component="R", sb_max_threshold=0.25):
    cov = float(record.INFO[f"{component}COVF"]) + float(
        record.INFO[f"{component}COVR"]
    )
    mincov = min(
        float(record.INFO[f"{component}COVF"]), float(record.INFO[f"{component}COVR"])
    )
    fq = mincov / cov
    minfq = min(1 - fq, fq)

    return alexSBFilter(cov, mincov, minfq, sb_max_threshold)


def get_vcf_strand_bias_filter(record, component="R"):
    cov = float(record.INFO[f"{component}COV"][0]) + float(
        record.INFO[f"{component}COV"][1]
    )
    mincov = min(
        float(record.INFO[f"{component}COV"][0]),
        float(record.INFO[f"{component}COV"][1]),
    )
    fq = mincov / cov
    minfq = min(1 - fq, fq)

    return alexSBFilter(cov, mincov, minfq)


# TODO: Move vals to cfg
# returns true if the record should be filtered
def alexSBFilter(cov, abs, fq, sb_max_threshold=0.25):
    if cov <= 10:
        return True
    elif cov <= 20:
        if abs < 5:
            return True
    elif cov <= 50:
        if abs < 10 and fq < sb_max_threshold:
            return True
    elif cov <= 100:
        if abs < 15 and fq < sb_max_threshold - 0.05:
            return True
    else:
        if fq < sb_max_threshold - 0.1:
            return True
    return False


def parseKmers(kmers, sequence, k):
    if len(sequence) < k:
        return
    for i in range(len(sequence) + 1 - k):
        kmer = str(sequence[i : i + k])
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1


def squashStrandedness(sign):
    if sign == ")":
        return "("
    else:
        return sign.upper()


# Parse pileup analysis file to a pythonesque structure
def parsePileup(pileupfile):
    pileupAnalysis = {}
    with open(pileupfile, "r") as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            pos = int(columns[0])
            ref = columns[1]
            spread = columns[2]
            spreadDict = {}
            for val in spread.split(";"):
                entry = val.split("=")
                spreadDict[squashStrandedness(entry[0])] = int(entry[1])
            pileupAnalysis[pos] = spreadDict
    return pileupAnalysis


# Parse pileup analysis file to a pythonesque structure
def parsePileupStrandAware(pileupfile):
    pileupAnalysis = {}
    with open(pileupfile, "r") as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            pos = int(columns[0])
            ref = columns[1]
            spread = columns[2]
            spreadDict = {}
            for val in spread.split(";"):
                entry = val.split("=")
                spreadDict[entry[0]] = int(entry[1])
            biases = []
            processed = set()
            for entry in spreadDict:
                mode = "low"
                if entry == ")" or entry.islower():
                    mode = "low"
                else:
                    mode = "high"

                if entry.upper() in processed:
                    continue
                else:
                    processed.add(entry.upper())

                for other in spreadDict:
                    if other == entry:
                        continue
                    if mode == "low" and other.lower() == entry:
                        biases += [abs(spreadDict[entry] - spreadDict[other])]
                        break
                    elif mode == "high" and other == entry.lower():
                        biases += [abs(spreadDict[entry] - spreadDict[other])]
                        break
            pileupAnalysis[pos] = (spreadDict, biases)
    return pileupAnalysis


def parsePileupStrandAwareLight(pileupfile):
    pileupAnalysis = {}
    with open(pileupfile, "r") as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            pos = int(columns[0])
            ref = columns[1]
            spread = columns[2]
            spreadDict = {}
            for val in spread.split(";"):
                entry = val.split("=")
                spreadDict[entry[0]] = int(entry[1])
            pileupAnalysis[pos] = spreadDict
    return pileupAnalysis


# returns the strand bias and the coverage for a given position and a pileup file
def parsePileupPosition(f, pos):
    with open(f, "r") as infile:
        for line in infile.read().splitlines():
            columns = line.split()
            if pos != columns[0]:  # column 0 contains the position
                continue  # we loop until we hit the position we are interested in
            spread = columns[2]  # column 2 contains the observed alleles
            spreadDict = {}
            for val in spread.split(";"):
                entry = val.split("=")
                spreadDict[entry[0]] = int(entry[1])
            bias = 0
            for entry in spreadDict:
                if entry.islower() or entry == ")":
                    bias -= spreadDict[entry]
                else:
                    bias += spreadDict[entry]

            return bias, sum(spreadDict.values())
        raise Exception("Can't find an entry for position {} in file {}".format(pos, f))


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(reversed([complement[n] for n in seq]))


def getStrandBias(pileupForPosition, altallele):
    altallele = "(" if altallele == "-" else altallele
    plus = 0
    total = 0
    for allele, count in pileupForPosition.items():
        if (allele.isupper() or allele == "(") and allele == altallele:
            plus += count
        if squashStrandedness(allele) == altallele:
            total += count
    return plus / total if total != 0 else 0


def getCoverage(pileupForPosition, altallele):
    altallele = "(" if altallele == "-" else altallele
    total = 0
    for allele, count in pileupForPosition.items():
        if squashStrandedness(allele) == altallele:
            total += count
    return total


def getTotalCoverage(pileupForPosition):
    return sum(pileupForPosition.values())


def getMinorStrandAbs(pileupForPosition, altallele):
    altallele = "(" if altallele == "-" else altallele
    total = 0
    for allele, count in pileupForPosition.items():
        if squashStrandedness(allele) == altallele and (
            allele == ")" or allele.islower()
        ):
            total += count
    return total


def getMinorStrandFrequency(pileupForPosition, altallele):
    sb = getStrandBias(pileupForPosition, altallele)
    return min(1 - sb, sb)


def getAlleleFrequency(pileupForPosition, altallele):
    altallele = "(" if altallele == "-" else altallele
    alleleCount = 0
    for allele, count in pileupForPosition.items():
        if squashStrandedness(allele) == altallele:
            alleleCount += count
    total = sum(pileupForPosition.values())
    return alleleCount / total
