import csv

from collections import defaultdict, Counter

from shared import *

def main(pangenome_path, reads_mapping, node2pos_path, rvt_threshold, th_het, output_path):

    # Get ref information
    node2pos = dict()
    node2ori = dict()
    ref_path = list()
    
    with open(node2pos_path) as fh:
        reader = csv.reader(fh, delimiter=",")
        next(reader)
        for row in reader:
            node2pos[row[0]] = int(row[1])
            node2ori[row[0]] = row[2]
            ref_path.append(row[0])
            
    ref_nodes = set(ref_path)
            
    # Get reads mappings
    node2cov = Counter()
    paths = defaultdict(list)
    parse_gaf(reads_mapping, paths, node2cov)

    # Found path didn't follow reference path
    variants = Counter() # [begin, variant path, end] -> list(reads)
    for path in paths:
        path_no_strand = [node[1:] for node in path]

        for node in path_no_strand:
            node2cov[node] += len(paths[path])
        
        if any(node in ref_nodes for node in path_no_strand): # Some node of path match with reference node
            if path_no_strand and not set(path_no_strand).issubset(ref_nodes): # Path isn't empty and some node isn't in reference path
                diff_pos = pos_of_diff(list(path_no_strand), list(ref_path))
                for diff_pos in pos_of_diff(list(path_no_strand), list(ref_path)):
                    if diff_pos is None:
                        continue
                    
                    variants[tuple(path[diff_pos[0]:diff_pos[1]])] += len(paths[path])

    with open("node2cov.csv", "w") as fh:
        print("node,coverage", file=fh)
        for (node, cov) in node2cov.items():
            print("{},{}".format(node, cov), file=fh)
        
    # set data required by vcf format
    ref_name = "MN908947.3"
    ref_length = 29903

    # get the sequence associate to node
    node2seq = get_node2seq(pangenome_path)

    final = dict()
    # write variant

    for (variant, count) in variants.items():
        variant_oris = [node[0] for node in variant]
        variant_nodes = [node[1:] for node in variant]

        ref_index = (ref_path.index(variant[0][1:]), ref_path.index(variant[-1][1:]))
        if ref_index[0] > ref_index[1]:
            ref_index = (ref_index[1], ref_index[0])
                
            variant_oris = ['>' if ori == '<' else '<' for ori in variant_oris]
            variant_nodes = list(reversed(variant_nodes))
                
        reference = ref_path[ref_index[0]:ref_index[1] + 1]

        reference_nodes = reference
        reference_oris = [node2ori[node] for node in reference_nodes]

        variant_seq = sequence_from_node(variant_oris, variant_nodes, node2seq, variant[0][0])
        reference_seq = sequence_from_node(reference_oris, reference_nodes, node2seq, variant[0][0])

        if variant_seq == reference_seq:
            continue
            
        pos = node2pos[reference_nodes[0]]
        ref_cov = min([node2cov[node] for node in reference_nodes])

        key = (variant_seq, reference_seq, pos, tuple(variant_nodes))
        if key in final:
            new_value = (final[key][0] + count, max(ref_cov, final[key][1]))
        else:
            new_value = (count, ref_cov)
        final[key] = new_value

    with open(output_path, "w") as fh:
        vcf_header(fh, ref_name, ref_length)

        for ((v_seq, r_seq, pos, nodes), (count, ref_cov)) in final.items():
            rvt = count / (count + ref_cov)

            if th_het <= rvt <= 1-th_het:
                v_seq = substituteAmbiguousBases(r_seq,v_seq)

            if rvt >= rvt_threshold:
                print("{}\t{}\t.\t{}\t{}\t.\t.\tVCOV={};RCOV={};RVT={};VARIANT_PATH={}".format(ref_name, pos + 1, r_seq, v_seq, count, ref_cov, rvt, ",".join(nodes)), file=fh)


#substitutes bases with their ambiguous counterparts if the threshold for heterozygosity is met for equal length non-SV vars
def substituteAmbiguousBases(r_seq,v_seq):
    if len(r_seq) == len(v_seq):
        ret = ''
        for a,b in zip(r_seq,v_seq):
            ambiguityChar = ambiguousBase(frozenset({a,b}))
            if ambiguityChar == None:
                return v_seq
            else:
                ret += ambiguityChar
        return ret
    else: #Do not substitute
        return v_seq
            
def vcf_header(fh, ref_name, length):
    print("##fileformat=VCFv4.2", file=fh)
    print("##INFO=<ID=VCOV,Number=1,Type=Integer,Description=\"Coverage of variant path\">", file=fh)
    print("##INFO=<ID=RCOV,Number=1,Type=Integer,Description=\"Coverage of reference path\">", file=fh)
    print("##INFO=<ID=RVT,Number=1,Type=Float,Description=\"Ratio between coverage of variant and total coverage\">", file=fh)
    print("##INFO=<ID=VARIANT_PATH,Number=1,Type=String,Description=\"Id of variant node in pangenome graph\">", file=fh)
    print("##contig=<ID={},length={}>".format(ref_name, length), file=fh)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=fh)

        
def sequence_from_node(oris, nodes, node2seq, previous_ori):
    seq = ""
    for ori, node in zip(oris, nodes):
        if ori == ">":
            seq += node2seq[node]
        else:
            seq += rev_comp(node2seq[node])
            
    return seq
                             
def pos_of_diff(path, reference):
    begin = path[0]
    
    try:
        index_ref = reference.index(begin)
    except ValueError:
        index_ref = None

    ref_set = set(reference)
    for i_not_ref in (index for (index, node) in enumerate(path) if node not in set(reference)):
        begin_break = None
        # found reference node before variant
        for i in range(i_not_ref, 0, -1):
            if path[i] in ref_set:
                begin_break = i
                break

        # found reference node after variant
        end_break = None
        for i in range(i_not_ref, len(path)):
            if path[i] in ref_set:
                end_break = i
                break

        if index_ref is None or begin_break is None or end_break is None:
            yield None
        else:
            yield (begin_break, end_break + 1)
    
if "snakemake" in locals():
    main(snakemake.input["pangenome"], snakemake.input["reads"], snakemake.input["node2pos"], float(snakemake.params["rvt"]),float(snakemake.params["th_het"]), snakemake.output["variant"])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]),float(sys.argv[5]), sys.argv[5])
