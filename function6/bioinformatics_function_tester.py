import pandas as pd
import numpy as np
from datetime import datetime
import os
import subprocess
from Bio import SeqIO

#By James C. Hu

def fq2fa(fq_path:str, fa_path:str) -> None:
    with open(fq_path, 'r') as f1, open(fa_path, 'w') as f2:
        sequences = SeqIO.parse(f1, 'fastq')
        count = SeqIO.write(sequences, f2, 'fasta')
        print(f'Converted {count} records')
    return None


fq2fa('/fastq/sample_S65_R1_001.fastq', '/fastq/sample_S65_R1_001.fasta')
