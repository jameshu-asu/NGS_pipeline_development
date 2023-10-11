from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt

#By James C. Hu

def tree_builder_fa(alined_fasta_file:str) -> None:
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
    tree = constructor.nj(distance_matrix) # options: nj or upgma
    print('\nPhylogenetic Tree\n==========================')
    Phylo.draw_ascii(tree)
    Phylo.write(tree, 'test_tree.nwk', 'newick')
    plt.rc('font', size=6)
    fig, ax = plt.subplots(figsize=(7, 9))
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('tree.png', dpi=200)
    plt.show()
    return None

Tree_infile = 'aligned.fasta'
tree_builder_fa(Tree_infile)











