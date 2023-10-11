from HAdV_iVar_helper import *
import time, os

#By James C. Hu

run = 365

ivar_directory_setup(run)
os.chdir('../')
run_samtools_coverage()
run_bedtools_coverage()
run_iVar_varaints('MN908947.3.fasta', 'MN908947.3.gff3')

os.chdir(f'iVar_coverage/bedtools_coverage')
with open(f'bedtools_coverage_summary.txt', 'w') as f:
    f.write(f'Seq_ID \t Total_nts \t Nts_with_reads \t Coverage \n ')
    for file in os.listdir():
        if file.endswith('bedtools_coverage.txt'):
            f.write(f'{file[:7]} \t {bedtools_depth_coverage_check(file, 1)[1]} \t {bedtools_depth_coverage_check(file, 1)[0]} \t {bedtools_depth_coverage_check(file, 1)[2]} \n')

os.chdir(f'../samtools_coverage')
final_samtools_df = pd.DataFrame()
for file in os.listdir():
    if file.startswith('I'):
        inum = file.split('_')[0]
        df = pd.read_csv(file, sep='\t')
        df['#rname'] = inum
        df = df.set_index('#rname')
        final_samtools_df = pd.concat([final_samtools_df, df])
final_samtools_df.to_csv('../../final_outputs/samtools_coverage_summary.txt')

os.chdir(f'../../iVar_tsv')
final_iVar_df = pd.DataFrame()
for file in os.listdir():
    if file.endswith('tsv'):
        df_out = iVar_variant_search(file)
        final_iVar_df = pd.concat([df_out, final_iVar_df])
final_iVar_df.to_csv('../final_outputs/Filtered_combined_iVar_output.csv')

os.chdir('..')
for file in os.listdir():
    if file.startswith('MN908947'):
        os.remove(file)
    if file.endswith(('_tsv', '_coverage', '_outputs')):
        shutil.move(file, f'iVar_deployable/{run}_iVar_analysis')


