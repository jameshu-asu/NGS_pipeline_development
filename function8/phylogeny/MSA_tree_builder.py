from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Align.Applications import MafftCommandline
from Bio.Phylo.Consensus import bootstrap_consensus
from Bio.Phylo.Consensus import majority_consensus
import matplotlib.pyplot as plt
import time
import sys


def MSA_mafft(infile: str) -> None:
    '''This funciton perfrom MSA on combined consensus file using MAFFT
    input: Combined consenssus fasta file.
    output: aligned.fasta
    '''
    mafft_exe = '/Users/chaoyihu/anaconda3/bin/mafft'
    mafft_cline = MafftCommandline(mafft_exe, input=infile)
    stdout, stderr = mafft_cline()
    with open('aligned.fasta', 'w') as handle:
        handle.write(stdout.upper())
    return None


def tree_builder_fa(alined_fasta_file: str) -> None:
    '''Builds Phylogenetic tree using alignment file generated during MSA
    input: Aligned fasta file
    output: Phylogenetic tree png, Newick format tree file.
    '''
    aligned_fasta = AlignIO.read(alined_fasta_file, 'fasta')
    print(aligned_fasta)
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(aligned_fasta)
    print('\nDistance Matrix\n==========================')
    print(distance_matrix)
    constructor = DistanceTreeConstructor()
    # options: nj (neighbor-joining) or upgma (Unweighted Pair Group Method with Arithemtic Mean)
    tree = constructor.nj(distance_matrix)
    print('\nPhylogenetic Tree\n==========================')
    Phylo.draw_ascii(tree)
    Phylo.write(tree, 'test_tree.nwk', 'newick')
    plt.rc('font', size=5)
    fig, ax = plt.subplots(figsize=(7, 9))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('tree.png', dpi=300)
    # plt.show()
    return None


def bootstrap_consensus_tree(aligned_fasta_file: str, replicates: int) -> None:
    '''Builds consensus tree based on n bootstrapped replicate trees.
    input: aligned fasta
    output: bootstrapped n replicate consensus tree.
    '''
    aligned_fasta = AlignIO.read(aligned_fasta_file, 'fasta')
    calculator = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor(calculator)
    tree = bootstrap_consensus(aligned_fasta, replicates, constructor, majority_consensus)
    Phylo.write(tree, 'consensus_tree.nwk', 'newick')
    plt.rc('font', size=5)
    fig, ax = plt.subplots(figsize=(7, 9))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('bootstrap_consensus_tree.png', dpi=300)
    Phylo.draw_ascii(tree)


MSA_infile = './consensus_fastas/combined_trimmedConsensus.fasta'
Tree_infile = 'aligned.fasta'


with open('timer_log.txt', 'w') as f:
    f.write(f'Process\tElapsed_time\n')

with open('timer_log.txt', 'a') as f:
    start = time.time()
    MSA_mafft(MSA_infile)
    end = time.time()
    delta = end - start
    f.write(f'MAFFT\t{delta}\n')

with open('timer_log.txt', 'a') as f:
    start = time.time()
    tree_builder_fa(Tree_infile)
    end = time.time()
    delta = end - start
    f.write(f'Tree_builder\t{delta}\n')

with open('timer_log.txt', 'a') as f:
    start = time.time()
    bootstrap_consensus_tree(Tree_infile, 100)
    end = time.time()
    delta = end - start
    f.write(f'Bootstrap_consensus_tree\t{delta}\n')
