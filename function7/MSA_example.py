from Bio.Align.Applications import MafftCommandline

#By James C. Hu

def MSA_mafft(infile: str) -> None:
    '''This function will:
    1) Perfrom MSA on combined consensus files.
    MSA: MAFFT
    
    input: Combined consenssus fasta file.
    output: aligned.fasta
    '''
    mafft_exe = '/Users/account/anaconda3/bin/mafft'
    mafft_cline = MafftCommandline(mafft_exe, input=infile)
    stdout = mafft_cline()
    with open ('aligned.fasta', 'w') as handle:
        handle.write(stdout.upper())
    return None

MSA_infile = './consensus_fastas/combined_trimmedConsensus.fasta'
MSA_mafft(MSA_infile)







