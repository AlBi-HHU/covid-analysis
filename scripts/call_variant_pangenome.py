log = open(snakemake.log[0], 'w')

graph_path = snakemake.input["graph"]
reads_path = snakemake.input["reads_mapping"]
reference_path = snakemake.input["reference_mapping"]

output_path = snakemake.output["variant"]

def parse_gaf(path, storage):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            tigs = re.findall(r"([<|>][^<>]+)", row[5])
        
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
        if path[i] == reference[index_ref + i]:
            end_break = i
            break

        
    if begin_break is None or end_break is None:
        return None
    else:
        begin_break = begin_break - 1
        end_break = end_break + 1
        return ((begin_break + index_ref, end_break + index_ref), (begin_break, end_break)) 

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))

def sequence_from_node(nodes, node2seq):
    sequence = node2seq[nodes[0][1:]]
    for node in nodes[1:]:
        seq = node2seq[node[1:]]

        if sequence[-20:] != seq[:20]:
            seq = rev_comp(seq)

        sequence += seq[20:]

    return sequence
        
import re
import csv

from collections import defaultdict

paths = defaultdict(list)
parse_gaf(reads_path, paths)

ref_paths = defaultdict(list)
parse_gaf(reference_path, ref_paths)
ref_path = next(iter(ref_paths.keys()))
ref_nodes = set(*ref_paths)

intresting_node = {node[1:] for node in ref_nodes}
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
                    paths[path]
                )

                intresting_node.update({v[1:] for v in variant[0]})
                intresting_node.update({v[1:] for v in variant[1]})
                
                variants.append(variant)

node2seq = dict()
with open(graph_path) as graph_fh:
    reader = csv.reader(graph_fh, delimiter='\t')
    for row in reader:
        if row[0] == "S" and row[1] in intresting_node:
            node2seq[row[1]] = row[2]

if not all(node[1:] in intresting_node for node in ref_path):
    assert False, "Error no set"
            
if not all(node[1:] in node2seq for node in ref_path):
    assert False, "Error no seq"
            
node2ref_pos = dict()
cumu_len = 0
for node in ref_path:
    node2ref_pos[node[1:]] = cumu_len
    cumu_len += len(node2seq[node[1:]])


with open(output_path, "w") as fh:
    print("##fileformat=VCFv4.2", file=fh)
    print("##INFO=<ID=CO,Number=1,Type=Integer,Description=\"Coverage of variant\">", file=fh)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", file=fh)

    for variant in variants:
        reference = sequence_from_node(variant[0], node2seq)
        read = sequence_from_node(variant[1], node2seq)
        pos = node2ref_pos[variant[0][0][1:]]
        coverage = len(variant[2])
        
        print("0\t{}\t.\t{}\t{}\t.\t.\tCO={}".format(pos, reference, read, coverage), file=fh)
    
