import pandas as pd
import numpy as np
import click as ck
from Bio import SeqIO
import os
from pathlib import Path
import gzip

ck.command()
def main():
    genomes = {
        'ksa001.1': '../t2t/data/ksa001.1v1.0.0.fa',
        'ksa001.2': '../t2t/data/ksa001.2v1.0.0.fa',
    }
    with open('data/input_ksa001.fa', 'w') as w:
        for key, seq_file in genomes.items():
            with open(seq_file, 'r') as f:
                seqs = SeqIO.parse(f, 'fasta')
                for rec in seqs:
                    rec.id = key + '|' + rec.id
                    SeqIO.write(rec, w, 'fasta')
            
if __name__ == '__main__':
    main()
