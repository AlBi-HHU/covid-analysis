from shared import get_node2seq
sys.path.append("scripts") #Hackfix but results in a more readable scripts folder structure
from shared import * #TODO: import only required modules
import csv
import re
from collections import defaultdict
import json

def main(alignment,pangenome,output):


    def get_node2seq(graph_path):
        node2seq = dict()

        with open(graph_path) as graph_fh:
            reader = csv.reader(graph_fh, delimiter='\t')
            for row in reader:
                if row[0] == "S" and row[1]:
                    node2seq[row[1]] = row[2]

        return node2seq

    matchesPerNode = defaultdict(int)

    # get the sequence associate to node
    node2seq = get_node2seq(pangenome)


    with open(alignment, 'r') as alignmentFile, open(snakemake.log[0],'w') as logFile, open(output,'w') as outFile:
        # For each line in the alignment file
        for line in alignmentFile.read().splitlines()[1:]:
            # We split based on whitespace
            data = line.split()

            # print('Alignment Block Length: {}'.format(int(data[10])))
            # print('Processing Read: {}'.format(data[0]))

            # Split Read, the first Symbol in path is only used as an indicator for the orientation
            path = list(data[5][1:].split('<')[::-1] if data[5][0] == '<' else data[5][1:].split('>'))

            # Not relevant: Print Path
            # print(path)

            pathLength = int(data[6])
            # print(int(data[7]),int(data[8]))

            pathStart = int(data[7]) if data[5][0] == '>' else pathLength - int(data[8])  # +1#SEE GAF SPECIFICATION
            pathEnd = int(data[8]) if data[5][0] == '>' else pathLength - int(data[7])

            # queryStart = int(data[2])
            # queryEnd = int(data[3])
            # print('Path Length: {}, Start on Path: [ {} End on Path: ) {}'.format(pathLength,pathStart,pathEnd))
            # print('Query Start: {}, Query End: {}'.format(queryStart,queryEnd))

            cigar = data[-1].split(':')[2]
            cigarSplit = re.findall('(\d+[MNDI])', cigar)
            if cigarSplit[-1][-1] == 'I':  # If the last instruction is an I, discard it, it has no information
                cigarSplit = cigarSplit[:-1]

            # Some Sanity Checking, remove this later TODO
            estimateTotal = 0
            for cigarInstruction in cigarSplit:
                if cigarInstruction[-1] != 'I':
                    estimateTotal += int(cigarInstruction[:-1])

            # print('Cigar covers {} bases'.format(estimateTotal))

            estimateTotalNodeLength = 0
            for node in path:
                estimateTotalNodeLength += len(node2seq[node])

            # print('Nodes on Path cover {} bases'.format(estimateTotalNodeLength))

            currentNode = path.pop(0)
            currentNodeLength = len(node2seq[currentNode])
            positionOnNode = 0  # local nucleotide position on the current node
            positionOnPath = 0  # global nucleotide position on the path

            ####### Just skip everything in the beginning noone cares

            # print('Processing {} bases before the alignment starts ... '.format(pathStart))
            for _ in range(pathStart):
                positionOnPath += 1
                positionOnNode += 1
                assert not (positionOnNode > currentNodeLength - 1)
                '''
                #went over the node boundary
                print('Finished node BEFORE cigar string started!{}'.format(currentNode))
                currentNode = path.pop(0)
                currentNodeLength = len(node2seq[currentNode])
                positionOnNode = 0
                '''

                ################ HERE STARTS CIGAR

            # print('Processing CIGAR String starting at path pos: {} (Node {} Processed: {})'.format(positionOnPath,currentNode,positionOnNode))

            # CigarInstruction = 13M
            for cigarInstructionIndex, cigarInstruction in enumerate(cigarSplit):

                # In this case 13
                cigarLength = int(cigarInstruction[:-1])

                # print('Processing CI {} (Instruction: {} with Length: {})'.format(cigarInstruction,cigarInstruction[-1],cigarLength))

                for cigarCounter in range(cigarLength):  # do x times where x is the cigar instruction count

                    # The info we want to have, how many Ms fit to a node
                    matchesPerNode[currentNode] += (cigarInstruction[-1] == 'M')

                    if cigarInstruction[-1] != 'I':  # don't move on insertion
                        positionOnPath += 1
                        positionOnNode += 1

                    # print('Now on position: {}'.format(positionOnPath))
                    # check if node is completed
                    if positionOnNode > currentNodeLength - 1:

                        if (cigarCounter == cigarLength - 1) and (cigarInstructionIndex == len(cigarSplit) - 1):
                            # print('Using Escape Condition')
                            pass
                        else:
                            '''
                            print('Moving to next node ... {}/{} of current CI: {} processed  '.format(
                                cigarCounter+1,cigarLength,cigarInstruction[-1])
                            )     
                            '''
                            currentNode = path.pop(0)
                            currentNodeLength = len(node2seq[currentNode])
                            positionOnNode = 0

            # print('Done with all CIGAR Stuff, GlobalPosition (post Alignment): {}'.format(positionOnPath))

            if cigarInstructionIndex != len(cigarSplit) - 1:
                logFile.write('Warning in Read: {}'.format(data[0]))
                logFile.write('Did not handle all cigars: {}/{}'.format(cigarInstructionIndex + 1, len(cigarSplit)))
            elif cigarLength < cigarCounter:
                logFile.write('Warning in Read: {}'.format(data[0]))
                logFile.write('Did not handle the final cigar comepletely: {}/{}'.format(cigarCounter, cigarLength))
            elif positionOnPath != pathEnd:
                logFile.write('Warning in Read: {}'.format(data[0]))
                logFile.write('Semi-Critical Error: Did not finish the path: {}/{}, currentNode: {}'.format(positionOnPath, pathEnd,
                                                                                                    currentNode))
            elif len(path) != 0:
                logFile.write('Warning in Read: {}'.format(data[0]))
                logFile.write('Critical Error: Path is not empty!!! : {}'.format(
                    ['{} (length:{})'.format(n, len(node2seq[currentNode])) for n in path]))
            else:
                pass
                # print('No problems processing this read!')

    json.dump(outFile,matchesPerNode)

if "snakemake" in locals():
    main(snakemake.input['alignment'], snakemake.input['pangenome'],snakemake.output[0])

else:
    import sys

    main(sys.argv[1], sys.argv[2],sys.argv[3])
