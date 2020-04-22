log = open(snakemake.log[0], 'w')

graph_path = snakemake.input['graph']
mapping_path = snakemake.input['mapping']

tig2cov_path = snakemake.output['tig2cov']
link2cov_path = snakemake.output['link2cov']

min_cov_tig = snakemake.params['min_cov_tig']
min_cov_link = snakemake.params['min_cov_link']

import re
import csv

from collections import defaultdict, Counter

tig2cov = defaultdict(list)
link2cov = defaultdict(list)

with open(mapping_path) as map_fh:
    reader = csv.reader(map_fh, delimiter='\t')
    for row in reader:
        tigs = re.split(r"[<|>]+", row[5])

        iterator = iter(tigs[1:]) # first element of list is empty
        prev = next(iterator)
        tig2cov[prev].append(row[0])

        for tig in iterator:
            tig2cov[tig].append(row[0])
            link2cov[frozenset((prev, tig))].append(row[0])
            
            prev = tig

            
with open(tig2cov_path, "w") as tig2cov_fh:
    for tig, cov in tig2cov.items():
        print("{},{}".format(tig, len(cov)), file=tig2cov_fh)
    
with open(link2cov_path, "w") as link2cov_fh:
    for link, cov in link2cov.items():
        link = sorted((link))
        print("{},{}".format(link, len(cov)), file=link2cov_fh)

        
print("histogram of tig coverage:", file=log)
tig_cov_histo = Counter(map(lambda x: len(x), tig2cov.values()))
for k in sorted(tig_cov_histo.keys()):
    print("{},{}".format(k, tig_cov_histo[k]), file=log)

print("histogram of link coverage:", file=log)
link_cov_histo = Counter(map(lambda x: len(x), link2cov.values()))
for k in sorted(tig_cov_histo.keys()):
    print("{},{}".format(k, tig_cov_histo[k]), file=log)
