
import re
import csv

degenerate = {
    frozenset(('A')): 'A',
    frozenset(('C')): 'C',
    frozenset(('G')): 'G',
    frozenset(('T')): 'T',
    frozenset(('A', 'G')): 'R',
    frozenset(('C', 'T')): 'Y',
    frozenset(('G', 'C')): 'S',
    frozenset(('A', 'T')): 'W',
    frozenset(('G', 'T')): 'K',
    frozenset(('A', 'C')): 'M',
    frozenset(('C', 'G', 'T')): 'B',
    frozenset(('A', 'G', 'T')): 'D',
    frozenset(('A', 'C', 'T')): 'H',
    frozenset(('A', 'C', 'G')): 'V',
    frozenset(('A', 'C', 'T', 'G')): 'N'
}

def ambiguousBase(frozenset):
    if frozenset in degenerate:
        return degenerate[frozenset]
    else:
        return None

def get_node2seq(graph_path):
    node2seq = dict()

    with open(graph_path) as graph_fh:
        reader = csv.reader(graph_fh, delimiter='\t')
        for row in reader:
            if row[0] == "S" and row[1]:
                node2seq[row[1]] = row[2]

    return node2seq

def parse_gaf(path, storage, node2base=None, edge2cov=None, node2seq=None):
    with open(path) as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            nodes = re.findall(r"([<|>][^<>]+)", row[5])
            begin_path = int(row[7])
            end_path = int(row[8])
            remain_base = end_path - begin_path

            if node2base is not None:
                first = True
                for node in nodes:
                    node_len = len(node2seq[node[1:]])
                    if first:
                        first = False
                        node2base[node[1:]] += node_len - begin_path
                        remain_base -= node_len
                    elif remain_base > node_len:
                        node2base[node[1:]] += node_len
                        remain_base -= node_len
                    else:
                        node2base[node[1:]] += remain_base
                        break

            if edge2cov is not None:
                for i in range(0, len(nodes) - 1):
                    edge2cov[frozenset((
                        nodes[i][1:],
                        nodes[i+1][1:]
                    ))] += 1

            storage[tuple(nodes)].append(row[0])


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def rev_comp(seq):
    return "".join(reversed([complement[n] for n in seq]))
