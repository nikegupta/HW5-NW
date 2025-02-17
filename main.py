# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    #inialize class
    nw = NeedlemanWunsch('substitution_matrices/BLOSUM62.mat',-10,-1)

    #dict to store alignments
    alignments = {'gallus_gallus': 0, 'mus_musculus': 0, 'balaeniceps_rex': 0, 'tursiops_truncatus': 0}
    keys = list(alignments.keys())

    #list of species to align to
    align_to = [gg_seq,mm_seq,br_seq,tt_seq]

    #run the alignments
    for i, seq in enumerate(align_to):
        nw.align(hs_seq,seq)
        alignments[keys[i]] = nw.alignment_score

    #print the species in list of highest alignment scores
    for i in sorted(alignments, key=alignments.get, reverse=True):
        print(i)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    for i in sorted(alignments, key=alignments.get, reverse=True):
        print(i + ' BRD2 alignment score with human BRD2: ' + str(alignments[i]))  

    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw.align(seq1,seq2)
    print(nw._align_matrix.tolist())

if __name__ == "__main__":
    main()
