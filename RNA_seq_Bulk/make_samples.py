'''
MAke samples.json file with sample and file names. https://raw.githubusercontent.com/slowkow/snakefiles/master/star_express/Snakefile
'''

import json
from glob import glob

fastqs = glob('/work/morrissy_lab/courtney/mGBM_modeling_clonality_RNA/Li*/*fastq.gz')
#fastqs = glob('/bulk/morrissy_bulk/INHOUSEDATA/Pericite_Biology/BaseCalls/PC-24CM-ECD2-12H*/*/*fastq.gz')
FILES= {}

SAMPLES = ["_".join(fastq.split('/')[-1].split('_',1)[:1]) for fastq in fastqs] 
#SAMPLES = ["_".join(fastq.split('/')[-1].split('_',3)[:3]) for fastq in fastqs]

#    FILES[mainsample] ={}
for sample in SAMPLES:
    # Change 'R1' and 'R2' to match the way your mate pairs are marked.
    mate1 = lambda fastq: sample in fastq and 'R1' in fastq
    mate2 = lambda fastq: sample in fastq and 'R2' in fastq
    FILES[sample] = {}
    FILES[sample]['R1'] = sorted(filter(mate1, fastqs))
    FILES[sample]['R2'] = sorted(filter(mate2, fastqs))

js = json.dumps(FILES, indent = 4, sort_keys=True)
open('samples_2.json', 'w').writelines(js)

