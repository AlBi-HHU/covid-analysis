log = open(snakemake.log[0], 'w')

threshold = snakemake.params['min_mapping_coverage']
path_prefix = snakemake.params['prefix']

graph_path = snakemake.input['graph']
mapping_path = snakemake.input['mapping']

output_path = snakemake.output['out_graph']

import re
import csv

from collections import defaultdict, Counter

tig2cov = defaultdict(list)
link2cov = defaultdict(list)

paths = defaultdict(list)

with open(mapping_path) as map_fh:
    reader = csv.reader(map_fh, delimiter='\t')
    for row in reader:
        tigs = re.split(r"(<|>)+", row[5])[1:]

        if len(tigs) <= 2:
            continue
        
        path = list()
        for (ori, tig) in zip(tigs[0::2], tigs[1::2]):
            if ori == '>':
                path.append(tig+'+')
            else:
                path.append(tig+'-')

        paths[tuple(path)].append(row[0])
        
with open(output_path, 'w') as out_fh:
    with open(graph_path, 'r') as in_fh:
        for line in in_fh:
            print(line, file=out_fh, end="")

    for i, (path, reads) in enumerate(paths.items()):
        if len(reads) > threshold:
            print("P\t{}{}\t{}\t*\tCO:i:{}\tRE:Z:{}".format(path_prefix, i, ",".join(path), len(reads), ",".join(reads)), file=out_fh)

        
