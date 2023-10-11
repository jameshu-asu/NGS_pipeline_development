import pandas as pd
import os
import subprocess
from pathlib import Path

#By James C. Hu

def file_transfer(input_file: pd.DataFrame) -> None:
    '''Input file MUST have the following columns:
    1) Seq_ID
    2) Instrument
    3) Instrument_run_number

    Note: Not sure if current pathing works for windows systems.

    This function transfers fastq files based on input criteria.
    1) Paths to sequencer directory using Instrument (sequencer name) and Instrument_run_number.
    2) Moves files that match Seq_IDs into respective working directory.
    '''
    sequencer_id = ('VL00126', 'VL00210')
    df = pd.read_csv(input_file).copy()
    df['Instrument'] = df['Instrument'].dropna()
    seq_id_list = list(set(df['Seq_ID'].to_list()))
    instrument_list = [value for value in df['Instrument'].to_list()
                       if str(value) != 'nan']
    instrument_run_number_list = [int(value) for value in df['Instrument_run_number'].to_list()
                                  if str(value) != 'nan']
    for instrument, instrument_run_number in zip(instrument_list, instrument_run_number_list):
        os.mkdir(f'{instrument}_{instrument_run_number}') # stitch to function1's directories
        os.mkdir(f'{instrument}_{instrument_run_number}/fastq')
        sequencer_output_directory = Path(f'/mnt/storage/data/{instrument}/outputDirectory/')
        run_directory_list = [directory for directory in os.listdir(sequencer_output_directory)
                              if any(match_criteria in directory for match_criteria in sequencer_id)
                              if directory.split('_')[2] == str(instrument_run_number)]
        for file in run_directory_list:
            fastq_directory_path = Path(f'/mnt/storage/data/{instrument}/outputDirectory/{file}/Analysis/1/Data/fastq')
            fastq_list = [fastq for fastq in os.listdir(f'/mnt/storage/data/{instrument}/outputDirectory/{file}/Analysis/1/Data/fastq')
                          if fastq.startswith('I')]
            for fastq in fastq_list:
                fastq_seq_id = fastq.split('_')[0]
                if fastq_seq_id in seq_id_list:
                    fastq_file_path = Path(f'/mnt/storage/data/{instrument}/outputDirectory/{file}/Analysis/1/Data/fastq/{fastq}')
                    subprocess.call([f'cp', f'{fastq_file_path}', f'{instrument}_{instrument_run_number}/fastq'])
    return None



file_transfer('input_file.csv')






