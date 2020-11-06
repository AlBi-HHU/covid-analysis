from Bio import SeqIO
from argparse import ArgumentParser
from itertools import combinations
from collections import defaultdict
import sys 

def compare_sequences(q, r, name):
    assert(len(q) == len(r))
    mms = []
    dels = []
    ins = []
    in_deletion = False
    in_insertion = False
    query = q
    ref = r
    query_ng = query.replace("-","")
    ref_ng = ref.replace("-","")

    ridx = 0
    qidx = 0
    sidx = 0
    for s1,s2 in zip(query,ref):
        if s1 == "-":
            sidx += 1
            if s2 != "-":
                ridx += 1
        else:
            break

    i = 0
    for s1,s2 in zip(query[::-1],ref[::-1]):
        if s1 == "-":
            i += 1
        else:
            break
    if i != 0:
        query = query[:-i]
        ref = ref[:-i]
    assert(len(query) == len(ref))


    for s1,s2 in zip(query[sidx:],ref[sidx:]):
        sidx+=1
        if s1 == "-":
            if s2 == "-":
                # no change in qidx and ridx
                if in_insertion:
                    ins.append((ridx,query_ng[qstart-1:qidx-1],ref_ng[ridx-1]))
                    in_insertion = False
                elif in_deletion:
                    dels.append((rstart,query_ng[qidx-1],ref_ng[rstart-1:ridx-1]))
                    in_deletion = False
            else: # s2 = [ACGT]
                ridx += 1
                if in_insertion:
                    ins.append((ridx-1,query_ng[qstart-1:qidx-1],ref_ng[ridx-2]))
                    in_insertion = False
                if not in_deletion:
                    in_deletion = True
                    rstart = ridx-1
        else: # s1 = [ACGT]
            qidx += 1
            if s2 == "-":
                if in_deletion:
                    dels.append((rstart,query_ng[qidx-2],ref_ng[rstart-1:ridx-1]))
                    #print(sidx)
                    in_deletion= False
                if not in_insertion:
                    in_insertion = True
                    qstart = qidx-1
                    
            else: # s2 =[ACGT]
                ridx += 1
                if in_insertion:
                    ins.append((ridx-1,query_ng[qstart-1:qidx-1],ref_ng[ridx-2]))
                    in_insertion = False
                elif in_deletion:
                    dels.append((rstart,query_ng[qidx-2],ref_ng[rstart-1:ridx-1]))
                    in_deletion= False
                if not(s1 == s2 or s1== 'N' or s2=='N'):
                    mms.append((ridx,s1,s2))
    return [mms, ins, dels]


def print_vcf(snps,ins, dels, vcffile, mincount=1):
    #snps = set(filter( lambda x: snps[x] >= mincount, snps))
    #ins = set(filter( lambda x: ins[x] >= mincount, ins))
    #dels = set(filter( lambda x: dels[x] >= mincount, dels))

    header="""\
##fileformat=VCFv4.1
##contig=<ID=1,length=29903>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT\
"""

    allvars = defaultdict(int)
    for x in snps:
        allvars[x] = snps[x]
    for x in ins:
        allvars[x] = ins[x]
    for x in dels:
        allvars[x] = dels[x]
    sortedvars = sorted(allvars, key= lambda x: x[0])
    #print(sortedvars)
    with open(vcffile, "w") as vcf:
        vcf.write(header + "\n")
        #for var,count in 
        for var in sortedvars:
            count = allvars[var]
            #CHROM POS      ID         REF   ALT    QUAL  FILTER   INFO                     FORMAT
            pos, alt, ref = var
            vcf.write("\t".join(["NC_045512.2", str(pos), ".", ref, alt, str(count),"PASS", ".", "GT", "1"]) + "\n")


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("msafile") 
    parser.add_argument("refid") 
    parser.add_argument("vcfout") 
    parser.add_argument("--blacklist") 
    args = parser.parse_args()
    ref = {}
    refseq = ""
    seqs = {}

    blacklisted_ids = set()
    if args.blacklist:
        with open(args.blacklist) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                blacklisted_ids.add(line.strip());

    verbose = False

    # Read Fasta Data
    for record in SeqIO.parse(args.msafile, "fasta"):
        if args.refid in record.id:
            ref[record.id] = str(record.seq).upper()
            refseq = str(record.seq).upper()
        elif record.id in blacklisted_ids:
            if verbose:
                print("filtered " + record.id)
            continue
        else:
            seqs[record.id] = str(record.seq).upper()



    tsnps = defaultdict(int)
    tins = defaultdict(int)
    tdels = defaultdict(int)
    number = 0
    for idx, seq in seqs.items():
        number += 1
        print(number/len(seqs),end="\r")
        snps, ins, dels = compare_sequences(seq, refseq,idx)
        #print("\t".join([idx, str(len(ins)), str(len(dels))]))
        if verbose:
            #print(snps)
            print(ins)
            print(dels)
            print("---")
        for snp in snps:
            tsnps[snp] +=1 
        for in1 in ins:
            tins[in1] +=1
        for del1 in dels:
            tdels[del1] +=1
    #print(tsnps)
    #print(tins)
    #print(tdels)
    print_vcf(tsnps,tins,tdels,args.vcfout)

