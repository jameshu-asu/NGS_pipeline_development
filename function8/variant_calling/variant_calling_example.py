import pandas as pd
import numpy as np
from typing import Tuple
import subprocess, os, shutil


def ivar_directory_setup(run) -> None:
    os.mkdir(f'{run}_iVar_analysis')
    shutil.copy('MN908947.3.fasta', '../')
    shutil.copy('MN908947.3.gff3', '../')
    os.chdir('../')
    subprocess.call(['mkdir', 'iVar_tsv', 'iVar_coverage', 'final_outputs'])
    os.chdir('iVar_coverage')
    subprocess.call(['mkdir', 'bedtools_coverage', 'samtools_coverage'])
    return None


def run_iVar_varaints(ref_gene_fa: str, reg_gene_gff3: str) -> None:
    for file in os.listdir():
        if file.endswith('.bam'):
            subprocess.call(f'samtools mpileup -d 0 -A -aa -Q 0 {file} | ivar variants -p {file[:7]}_iVAR -q 20 -t 0 -r {ref_gene_fa} -g {reg_gene_gff3}', shell=True)
    for file in os.listdir():
        if file.endswith('.tsv'):
            shutil.move(file, 'iVar_tsv')
    return None


def run_bedtools_coverage() -> None:
    for file in os.listdir():
        if file.endswith('.bam'):
            subprocess.call(f'genomeCoverageBed -d -ibam {file} > {file[:7]}_bedtools_coverage.txt', shell=True)
    for file in os.listdir():
        if file.endswith('_bedtools_coverage.txt'):
            shutil.move(file, 'iVar_coverage/bedtools_coverage')
    return None


def run_samtools_coverage() -> None:
    for file in os.listdir():
        if file.endswith('.bam'):
            subprocess.call(['samtools', 'coverage', f'{file}', '-o', f'{file[:7]}_samtools_coverage.txt'])
    for file in os.listdir():
        if file.endswith('_samtools_coverage.txt'):
            shutil.move(file, 'iVar_coverage/samtools_coverage')
    return None


def bedtools_avg_depth_at_loci(file) -> Tuple[float, float, float]:
    df = pd.read_csv(file, sep='\t', names=['REGION', 'NT_POS', 'READS' ])
    df_E180 = df[df['NT_POS'].between(22100, 22102, inclusive='both')].set_index('REGION')
    E180_avg_reads = df_E180['READS'].mean()
    df_K478 = df[df['NT_POS'].between(22994, 22996, inclusive='both')].set_index('REGION')
    K478_avg_reads = df_K478['READS'].mean()
    df_F486 = df[df['NT_POS'].between(23018, 23020, inclusive='both')].set_index('REGION')
    F486_avg_reads = df_F486['READS'].mean()
    return E180_avg_reads, K478_avg_reads, F486_avg_reads;


def bedtools_depth_coverage_check(file, threshold: int) -> Tuple[int, int, float]:
    df = pd.read_csv(file, sep='\t', names=['REGION', 'NT_POS', 'READS' ])
    nts_with_reads = (df['READS'] > threshold).sum()
    total_nts = len(df['READS'])
    coverage = nts_with_reads / total_nts
    return nts_with_reads, total_nts, coverage;


def iVar_variant_search(file) -> pd.DataFrame:
    df = pd.read_csv(file, sep='\t')
    df['REGION'] = file[:7]
    df_E180 = df[df['POS'].between(22100, 22102, inclusive='both')].set_index('REGION')
    df_K478 = df[df['POS'].between(22994, 22996, inclusive='both')].set_index('REGION')
    df_F486 = df[df['POS'].between(23018, 23020, inclusive='both')].set_index('REGION')
    df_out = pd.concat([df_E180, df_K478, df_F486])
    return df_out




