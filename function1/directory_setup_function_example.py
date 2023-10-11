import pandas as pd
import os
from itertools import zip_longest

#By James C. Hu
# Function to be imported from production script.
def generate_project_directory(input_file: pd.DataFrame) -> None:
    '''
    This script will generate the directory structure for new projects.

    Example structure

    Root_directory
    |- new_project_name
    |  |- report_files
    |  |  |- coverage_files
    |  |  |- alignment_files
    |  |  |- consensus_files
    |  |- sample_files
    |  |  |- fastq_files
    |  |  |- fasta_files
    |  |- log_files
    '''
    df = pd.read_csv(input_file, index_col=0).copy()

    directory_ids = df.index.to_list()

    parent = df.at['parent', 'Directory_name']
    os.mkdir(parent)

    main_branches_list = list(set([directory.split('_')[0] for directory in directory_ids if directory.startswith('b')]))
    main_branches_list.sort()

    for branch in main_branches_list:
        branch = df.at[branch, 'Directory_name']
        os.mkdir(f'{parent}/{branch}')


    sub_branches_list = list(set([directory for directory in directory_ids if '_' in directory]))
    sub_branches_list.sort()

    for main_branch_ID in main_branches_list:
        main_branch_name = df.at[main_branch_ID, 'Directory_name']
        print(main_branch_name)
        for sub_branch_ID in sub_branches_list:
            sub_branch_name = df.at[sub_branch_ID, 'Directory_name']
            if sub_branch_ID.split('_')[0] == main_branch_ID:
                print(f'{parent}/{main_branch_name}/{sub_branch_name}')
                os.mkdir(f'{parent}/{main_branch_name}/{sub_branch_name}')
    return None

# Funciton calling method in project script
generate_project_directory('input_file.csv')





