# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        #make arrays
        self._align_matrix = np.zeros((len(self._seqA) + 1, len(self._seqB) + 1)) 
        self._gapA_matrix = np.zeros((len(self._seqA) + 1, len(self._seqB) + 1))
        self._gapB_matrix = np.zeros((len(self._seqA) + 1, len(self._seqB) + 1))

        #fill first row
        for i in range(1,self._align_matrix.shape[1]):
            self._align_matrix[0,i] = self.gap_open + self.gap_extend * i
            self._gapA_matrix[0,i] = self.gap_open + self.gap_extend * i
            self._gapB_matrix[0,i] = self.gap_open + self.gap_extend * i

        #fill first column
        for i in range(1,self._align_matrix.shape[0]):
            self._align_matrix[i,0] = self.gap_open + self.gap_extend * i
            self._gapA_matrix[i,0] = self.gap_open + self.gap_extend * i
            self._gapB_matrix[i,0] = self.gap_open + self.gap_extend * i

        # TODO: Implement global alignment here
        for a_i in range(1,self._align_matrix.shape[0]): #iterate through rows, corresponds to seqA, start at 1 because 0th row/column is 0
           for b_i in range(1,self._align_matrix.shape[1]): #iterate through columns, corresponds to seqB
               
                #calculate the match score
                match_score = self._align_matrix[a_i-1,b_i-1] + self.sub_dict[(self._seqA[a_i-1],self._seqB[b_i-1])]
                
                #calculate gap scores
                self._gapA_matrix[a_i,b_i] = max(self._align_matrix[a_i-1,b_i] + self.gap_open+self.gap_extend, self._gapA_matrix[a_i-1,b_i] + self.gap_extend)
                self._gapB_matrix[a_i,b_i] = max(self._align_matrix[a_i,b_i-1] + self.gap_open+self.gap_extend, self._gapB_matrix[a_i,b_i-1] + self.gap_extend)

                #find the highest score
                self._align_matrix[a_i,b_i] = max(match_score, self._gapA_matrix[a_i,b_i], self._gapB_matrix[a_i,b_i])

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        #get alignemtn score from bottom corner
        self.alignment_score = self._align_matrix[len(self._seqA)][len(self._seqB)]

        #raw alignments of A/B (are backwards)
        align_A = []
        align_B = []

        #initialize indices to bottom corner
        a_i = self._align_matrix.shape[0] - 1
        b_i = self._align_matrix.shape[1] - 1

        while a_i or b_i > 0:

            #step diagonally if score matches blossum score
            if self._align_matrix[a_i][b_i] == self._align_matrix[a_i-1][b_i-1] + self.sub_dict[(self._seqA[a_i-1],self._seqB[b_i-1])]:
                align_A.append(self._seqA[a_i-1])
                align_B.append(self._seqB[b_i-1])
                a_i -= 1
                b_i -= 1
            
            #step up if matches gap A matrix
            elif self._align_matrix[a_i][b_i] == self._gapA_matrix[a_i][b_i]:
                align_A.append(self._seqA[a_i-1])
                align_B.append('-')
                a_i -= 1
            
            #step left if matches gap B matrix
            elif self._align_matrix[a_i][b_i] == self._gapB_matrix[a_i][b_i]:
                align_A.append('-')
                align_B.append(self._seqB[b_i-1])
                b_i -= 1

            #condition should not be possible but is here just in case
            else:
                raise ValueError('Backtracing failed')

        #reverse and collapse into strings the alignments
        self.seqA_align = ''.join(align_A[::-1])
        self.seqB_align = ''.join(align_B[::-1])

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
