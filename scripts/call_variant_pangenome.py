log = open(snakemake.log[0], 'w')

import re
import csv

from collections import defaultdict

from shared import rev_comp

def main(graph_path, reads_path, reference_path, node_pos_on_ref, output_path, ovl_len, rvt_threshold):

    node2cov = defaultdict(int)
    
    # Parse mapping
    paths = defaultdict(list)
    parse_gaf(reads_path, paths, node2cov)

    ref_paths = defaultdict(list)
    parse_gaf(reference_path, ref_paths, defaultdict(int))
    ref_path = next(iter(ref_paths.keys()))
    ref_nodes = set(*ref_paths)

    # Get sequence of each node
    node2seq = dict()
    with open(graph_path) as graph_fh:
        reader = csv.reader(graph_fh, delimiter='\t')
        for row in reader:
            if row[0] == "S" and row[1]:
                node2seq[row[1]] = row[2]

    # Get position of node on reference
    node2ref_pos = dict()
    cumu_len = 0
    for node in ref_path:
        node2ref_pos[node[1:]] = cumu_len
        cumu_len += len(node2seq[node[1:]]) - ovl_len

    with open(node_pos_on_ref, "w") as csv_fh:
        print("node,pos", file=csv_fh)
        for node, pos in sorted(node2ref_pos.items(), key=lambda x: x[1]):
            print("{},{}".format(node, pos), file=csv_fh)

    
    # Found variant path
    variants = list()
    for path in paths.keys():
        if any(node in ref_nodes for node in path):
            if not set(path).issubset(ref_nodes):
                if path is None:
                    continue
            
                diff_pos = pos_of_diff(list(path), list(ref_path))
                if diff_pos is not None:
                    variant = (
                        ref_path[diff_pos[0][0]:diff_pos[0][1]],
                        path[diff_pos[1][0]:diff_pos[1][1]],
                    )

                    variants.append(variant)                    

                    
    ref_name = next(iter(ref_paths.values()))[0]
    
    # Concatenate variants on same node
    variants.sort()
    simplify_variants = list()
    if len(variants) == 0:
        with open(output_path, "w") as fh:
            vcf_header(fh, ref_name)
        return
    
    pred = variants[0]

    for variant in variants[1:]:
        if (pred[0], pred[1]) == (variant[0], variant[1]):
            pass
        else:
            simplify_variants.append((pred[0], pred[1]))
           
            pred = variant

    simplify_variants.append((pred[0], pred[1]))

    # Write variant
    with open(output_path, "w") as fh:
        vcf_header(fh, ref_name)

        for variant in simplify_variants:
            reference = sequence_from_node(variant[0], node2seq, ovl_len)
            read = sequence_from_node(variant[1], node2seq, ovl_len)
            
            pos = node2ref_pos[variant[0][0][1:]]

            var_cov = min([node2cov[node] for node in variant[1]])
            ref_cov = min([node2cov[node] for node in variant[0]])

            rvt = var_cov / (var_cov + ref_cov)
            
            if rvt >= rvt_threshold:
                print("{}\t{}\t.\t{}\t{}\t.\t.\tVCOV={};RCOV={};RVT={}".format(ref_name, pos + 1, reference, read, var_cov, ref_cov, rvt), file=fh)


def vcf_header(fh, ref_name):
    print("##fileformat=VCFv4.2", file=fh)
    print("##INFO=<ID=VCOV,Number=1,Type=Integer,Description=\"Coverage of variant path\">", file=fh)
    print("##INFO=<ID=RCOV,Number=1,Type=Integer,Description=\"Coverage of reference path\">", file=fh)
    print("##INFO=<ID=RVT,Number=1,Type=Float,Description=\"Ratio between coverage of variant and total coverage\">", file=fh)
    print("##contig=<ID={},length=29903>".format(ref_name), file=fh)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=fh)

def parse_gaf(path, storage, node2cov):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            tigs = re.findall(r"([<|>][^<>]+)", row[5])

            for tig in tigs:
                node2cov[tig] += 1
            
            storage[tuple(tigs)].append(row[0])

def pos_of_diff(path, reference):
    begin = path[0]
    try:
        index_ref = reference.index(begin)
    except ValueError:
        return None

    begin_break = None
    for i in range(0, len(path)):
        if path[i] != reference[index_ref + i]:
            begin_break = i
            break

    end_break = None
    for i in range(i, len(path)):
        if index_ref + i >= len(reference):
            end_break = None
            break
        
        if path[i] == reference[index_ref + i]:
            end_break = i
            break

        
    if begin_break is None or end_break is None:
        return None
    else:
        begin_break = begin_break - 1
        end_break = end_break + 1
        return ((begin_break + index_ref, end_break + index_ref), (begin_break, end_break)) 

    
def sequence_from_node(nodes, node2seq, ovl_len):
    if nodes[0][0] == '>':
        sequence = node2seq[nodes[0][1:]]
    else:
        sequence = rev_comp(node2seq[nodes[0][1:]])
    for node in nodes[1:]:
        seq = node2seq[node[1:]]

        if sequence[-1 * ovl_len:] != seq[:ovl_len]:
            seq = rev_comp(seq)

        sequence += seq[ovl_len:]

    return sequence
        
main(snakemake.input["graph"], snakemake.input["reads_mapping"], snakemake.input["reference_mapping"], snakemake.output["node_pos_on_ref"], snakemake.output["variant"], snakemake.params["ksize"] - 1, snakemake.params["rvt_threshold"])
