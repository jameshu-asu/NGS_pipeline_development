import pandas as pd
import numpy as np
from datetime import datetime
import os
import subprocess
from Bio import SeqIO

#By James C. Hu
#______General_____#

def concat_files(input1_path:str, input2_path:str, output_name_path:str) -> None:
    '''Input: Any 2 files to be concatinated.
    COVIDSeq example: Merging of deduplicated merged and unmerged fastq files.

    Validated 7.25.23
    '''
    with open(f'{output_name_path}', 'w') as f:
        subprocess.Popen(['cat', f'{input1_path}', f'{input2_path}'], stdout=f)
    return None


def copy_file(target_path:str, destination_path:str) -> None:
    '''This function is a python wrapper for CLI cp command.

    Validated 7.25.23
    '''
    subprocess.call(['cp', f'{target}', f'{destination}'])
    return None


def fq_2_fa(fq_path:str, fa_path:str) -> None:
    '''Input: Fastq file path, fasta file path (inclide file name)
    Output: Fasta

    Validated 7.25.23
    '''
    with open(fq_path, 'r') as f1, open(fa_path, 'w') as f2:
        sequences = SeqIO.parse(f1, 'fastq')
        count = SeqIO.write(sequences, f2, 'fasta')
    return None


#_____bbtools_____#

def bbduk1(seq_id: str, script_path:str, r1_in:str, r2_in:str, ref:str, r1_out:str,
          r2_out:str, k_val:int, hdist_val:int, ktrim_val:str, qtrim_val:str,
          mink_val:int, trimq_val:int, min_len_val:int, min_avg_qual_val:int, remove_if_val:str) -> None:
    '''Input: Raw reads from paired-end fastq files.
    Output: Indexes/adapter trimmed raw reads fastq files.

    Validated 7.25.23
    '''
    with open(f'{seq_id}_adaptTrimQC_log.txt', 'w') as f:
        subprocess.call([f'{script_path}',
                         f'in1={r1_in}',
                         f'in2={r2_in}',
                         f'ref={ref}',
                         f'out1={r1_out}',
                         f'out2={r2_out}',
                         f'k={k_val}',
                         f'hdist={hdist_val}',
                         f'ktrim={ktrim_val}',
                         f'qtrim={qtrim_val}',
                         f'mink={mink_val}',
                         f'trimq={trimq_val}',
                         f'minlength={min_len_val}',
                         f'minavgquality={min_avg_qual_val}',
                         f'removeifeitherbad={remove_if_val}'], stderr=f)
    return None



def bbduk_qc_PhiX(bbduk_in1:str, bbduk_in2:str, bbduk_ref_path:str, bbduk_r1_out_path:str, bbduk_r2_out_path:str, bbduk_k=None, bbduk_hdist=None, bbduk_overwrite=None) -> None:
    '''Step2 is only for COVIDSeq when using PhiX as positive control.
    This qc step removes PhiX reads from the sample.
    '''
    subprocess.call(['bbduk.sh',
                f'in={bbduk_in1}',
                f'in2={bbduk_in2}',
                f'ref={bbduk_ref_path}',
                f'out={bbduk_r1_out_path}',
                f'out2={bbduk_r2_out_path}',
                f'k={bbduk_k}',
                f'hdist={bbduk_hdist}',
                f'overwrite={bbduk_overwrite}'])
    return None


def bbduk_dedupe(dedupe_Xmx=None, dedupe_in1:str, dedupe_in2:str, dedupe_out1_path:str, dedupe_out2_path:str, dedupe_csf=None, dedupe_overwrite=None, dedupe_minidentity:int) -> None:
    '''Input: Trimmed fastq
    Output: Deduplicated trimmed fastq
    '''
    subprocess.call(['dedupe.sh',
                f'-Xmx{dedupe_Xmx}',
                f'in1={dedupe_in1}'
                f'in2={dedupe_in2}'
                f'out={dedupe_out1_path}'
                f'outd={dedupe_out2_path}'
                f'csf={dedupe_csf}'
                f'overwrite={dedupe_overwrite}',
                f'minidentity={dedupe_minidentity}'])
    return None


def bbmerge(bbmerge_in:str, bbmerge_out1_path:str, bbmerge_out2_path:str)-> None:
    '''Input: Unmerged fastq reads
    Output: Merged overlapping paired-read fastq.
    '''
    subprocess.call(['bbmerge.sh',
                 f'in={bbmerge_in}'
                 f'out={bbmerge_out1_path}'
                 f'outu={bbmerge_out2_path}'])
    return None


#_____bowties2______#

def bowtie2(bowtie2_x_path:str, bowtie2_f_path:str, bowtie2_out_path:str) -> None:
    '''Local alignemnt to reference sequnce using bowtie2
    bowtie2_x_path: Path to reference sequences .fasta file.
    bowtie2_f_path: Input file to be mapped as qc'd, deduped, merged sample .fasta file.
    bowtie2_out_path: Mapped output as .sam file.
    '''
    subprocess.call(['bowtie2', '-x', f'{bowtie2_x_path}', '-f', f'{bowtie2_f_path}', '>', f'{bowtie2_out_path}'])
    return None


#_____samtools______#

def samtools_view(samtools_F:hex, samtools_in_path:str, samtools_out_path:str) -> None:
    '''Input: Any file that has reads that need to be extracted.
    Output: File with which the extracted reads will be combined.

    FLAGS:
    0x1         PAIRED          paired-end (or multiple-segment) sequencing technology
    0x2         PROPER_PAIR     each segment properly aligned according to the aligner
    0x4         UNMAP           segment unmapped
    0x8         MUNMAP          next segment in the template unmapped
    0x10        REVERSE         SEQ is reverse complemented
    0x20        MREVERSE        SEQ of the next segment in the template is reverse complemented
    0x40        READ1           the first segment in the template
    0x80        READ2           the last segment in the template
    0x100       SECONDARY       secondary alignment
    0x200       QCFAIL          not passing quality controls
    0x400       DUP             PCR or optical duplicate
    0x800       SUPPLEMENTARY   supplementary alignment
    '''
    subprocess.call(['samtools', 'view', '-h', '-F', f'{samtools_F}', f'{samtools_in_path}', '>', f'{samtools_out_path}'])
    return None


def samtools_sort(samtools_in_path:str, samtools_out_path:str) -> None:
    '''Input: SAM/BAM/CRAM files
    Output: Sorted SAM/BAM/CRAM files by feature (leftmost coordinate, read name, tag content, and mimimiser-based collation order).
    '''
    subprocess.call(['samtools', 'sort',  f'{samtools_in_path}', '>',  f'{samtools_out_path}'])
    return None


def samtools_index(samtools_in_path:str) -> None:
    '''Input: SAM/BAM/CRAM file
    Output: Indexed SAM/BAM/CRAM file
    '''
    subprocess.call(['samtools', 'index',  f'{samtools_in_path}'])
    return None


def samtools_idxstats(samtools_in_path:str, samtools_in_path:str) -> None:
    '''
    '''
    subprocess.call(['samtools', 'idxstats',  f'{samtools_in_path}', '>',  f'{samtools_out_path}'])
    return None







