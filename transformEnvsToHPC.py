#Replaces the channels with channels that work on HPC

import os

directory = 'envs'


for filename in os.listdir(directory):
    if filename.endswith('.yaml'):
        print('Processing file: {}'.format(filename))
        out = ''
        with open(os.path.join(directory, filename),'r') as filehandler:
                out = filehandler.read()
        out = out.replace('defaults','http://conda.repo.test.hhu.de/main')        
        out = out.replace('bioconda','http://conda.repo.test.hhu.de/bioconda')        
        out = out.replace('conda-forge','http://conda.repo.test.hhu.de/conda-forge')   
        out = out.replace('channels:','channels:\n  - nodefaults')        
        with open(os.path.join(directory, filename),'w') as filehandler:
                filehandler.write(out)  
    else:
        continue
