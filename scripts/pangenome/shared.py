
import re
import csv

def get_node2seq(graph_path):
    node2seq = dict()
    
    with open(graph_path) as graph_fh:
        reader = csv.reader(graph_fh, delimiter='\t')
        for row in reader:
            if row[0] == "S" and row[1]:
                node2seq[row[1]] = row[2]

    return node2seq

def parse_gaf(path, storage, node2cov=None):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            tigs = re.findall(r"([<|>][^<>]+)", row[5])
                    
            if node2cov:
                for tig in tigs:
                    node2cov[tig] += 1
            
            storage[tuple(tigs)].append(row[0])


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))
