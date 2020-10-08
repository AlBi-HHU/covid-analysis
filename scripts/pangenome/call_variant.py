import csv

from collections import defaultdict, Counter

from shared import parse_gaf, get_node2seq, rev_comp

def main(pangenome_path, reads_mapping, node2pos_path, rvt_threshold, output_path):

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
                if diff_pos:
                    variants[tuple(path[diff_pos[0]:diff_pos[1]])] += len(paths[path])

    with open("node2cov.csv", "w") as fh:
        print("node,coverage", file=fh)
        for (node, cov) in node2cov.items():
            print("{},{}".format(node, cov), file=fh)

    with open("clean_graph.gfa", "w") as ofh:
        with open(pangenome_path) as ifh:
            for line in ifh:
                if line.startswith("S"):
                    row = line.split("\t")
                    if row[1] in node2cov:
                        print(line, file=ofh, end="")
                elif line.startswith("L"):
                    row = line.split("\t")
                    if row[1] in node2cov and row[3] in node2cov:
                        print(line, file=ofh, end="")
                else:
                    print(line, file=ofh)
            
    # set data required by vcf format
    ref_name = "MN908947.3"
    ref_length = 29903

    # get the sequence associate to node
    node2seq = get_node2seq(pangenome_path)

    # write variant
    with open(output_path, "w") as fh:
        vcf_header(fh, ref_name, ref_length)

        for (variant, count) in variants.items():
            variant_oris = [node[0] for node in variant[1:-1]]
            variant_nodes = [node[1:] for node in variant[1:-1]]

            #print(variant)
            reference = ref_path[ref_path.index(variant[0][1:])+1:ref_path.index(variant[-1][1:])]
            #print("reference ", reference)

            reference_nodes = reference
            reference_oris = [node2ori[node] for node in reference_nodes]

            
            #print("variant ori", variant_oris)
            #print("variant node", variant_nodes)

            #print("reference ori", reference_oris)
            #print("reference node",reference_nodes)

            variant_seq = sequence_from_node(variant_oris, variant_nodes, node2seq, variant[0][0])
            reference_seq = sequence_from_node(reference_oris, reference_nodes, node2seq, variant[0][0])

            pos = node2pos[reference_nodes[0]]
            ref_cov = min([node2cov[node] for node in reference_nodes])

            rvt = count / (count + ref_cov)

            if rvt >= rvt_threshold:
                print("{}\t{}\t.\t{}\t{}\t.\t.\tVCOV={};RCOV={};RVT={};VARIANT_PATH={}".format(ref_name, pos + 1, reference_seq, variant_seq, count, ref_cov, rvt, ",".join(variant_nodes)), file=fh)
            
            #print(variant_seq)
            #print(reference_seq)
            #print(node2pos[reference_nodes[0]], len(reference_seq))

            
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
        return None

    begin_break = None
    for i in range(0, len(path)):

        if index_ref + i >= len(reference):
            continue
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

        return (begin_break, end_break)

if "snakemake" in locals():
    main(snakemake.input["pangenome"], snakemake.input["reads"], snakemake.input["node2pos"], float(snakemake.params["rvt"]), snakemake.output["variant"])
else:
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
