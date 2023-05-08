"""
    This python file creates the GENE_IDS required for 
    kallisto to run. 

    It uses tx2gene to read the GTF Files

    Edwin Huang
    edh021@cloud.ucsd.edu
"""

import argparse
import pandas as pd
import numpy as np
from gtfparse import read_gtf

def tx2gene(gtf_file):
    """
    Reads the gtf file, does transformations, and outputs 
    a t2g.csv for use in transcript2genes.R

    """

    df = read_gtf(gtf_file)
    print(f'Shape of the GTF file: {df.shape}')
    t2g = df[df['feature']=='transcript'][['transcript_id','gene_id','gene_name']]
    t2g['target_id'] = t2g['transcript_id']
    uniq = len(np.unique(t2g['transcript_id']))
    print(f'Number of unique transcript IDs: {uniq}')
    uniq_gene = len(np.unique(t2g['gene_id']))
    print(f'Number of unique gene IDs: {uniq_gene}')
    uniq_name = len(np.unique(t2g['gene_name']))
    print(f'Number of unique gene names: {uniq_name}')

    ## writing t2g to csv:
    if "v43" in gtf_file: ## using a human annotation
        t2g.to_csv('HUMAN_gencode.v43_t2g.csv',index=False)
    elif "vM32" in gtf_file: ## this is a mouse annotation
        t2g.to_csv('MOUSE_gencode.vM32_t2g.csv',index=False)
    
    return 'Finished'



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-gtf', '--transcript_gtf', help = 'The Transcript GTF file to use')
    args = parser.parse_args()

    species_type = args.transcript_gtf
    tx2gene(species_type)
    