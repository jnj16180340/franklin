#!/usr/bin/python3

import string
import itertools
import sys
import argparse
# import re

def dna_to_rna(sequence):
    '''Replaces "T" with "U" in string'''
    replacement_table = str.maketrans('ATCG','AUCG')
    return sequence.translate(replacement_table)

def rna_to_dna(sequence):
    '''Replaces "U" with "T" in string '''
    replacement_table = str.maketrans('AUCG','ATCG')
    return sequence.translate(replacement_table)

def reverse_complement(sequence):
    '''Generates reverse complement of a DNA sequence string '''
    replacement_table = str.maketrans('ATCG','TAGC')
    return sequence[::-1].translate(replacement_table)

# TODO: for testing
#assert reverse_complement('ATCG') == 'TAGC'[::-1], 'bad'

def print_orfs(seq,codon_table):
    ''' Finds open reading frames in DNA string, including nested ORFs, and prints them, along with RNA and protein translations'''
    # TODO: validate against BioPython, NCBI tools
    # TODO: return positions instead of ORF strings
    # TODO: This doesn't include the stop codons, which is technically wrong. To include stop codons, implement a buffered iterator so that itertools.takewhile returns the item on which it stops.
    stop_codons = [rna_to_dna(k) for k, v in codon_table.items() if v == '*']
    start_codons = [rna_to_dna(k) for k, v in codon_table.items() if v == 'M']
    orfs = []
    # if we care about positions:
    #frame1 = zip(itertools.count(),(''.join(k) for k in zip(seq[0::3],seq[1::3],seq[2::3])))
    #frame2 = zip(itertools.count(),(''.join(k) for k in zip(seq[1::3],seq[2::3],seq[3::3])))
    #frame3 = zip(itertools.count(),(''.join(k) for k in zip(seq[2::3],seq[3::3],seq[4::3])))
    
    def chunk3frames(frnum): 
        '''Split up DNA sequence string into triplets, offset by frame '''
        return (''.join(k) for k in zip(seq[0+frnum::3],seq[1+frnum::3],seq[2+frnum::3]))
    
    for frame in map(chunk3frames,range(3)):
        exhausted = False # Are there no more ORFs to find?
        passthrough = itertools.dropwhile(lambda l: l not in start_codons, frame)
        while exhausted is False:
            passthrough, process = itertools.tee(passthrough)
            result = itertools.takewhile(lambda l: l not in stop_codons, process) # this omits the stop codon
            new_orf = list(result)
            passthrough = itertools.dropwhile(lambda l: l not in start_codons, itertools.islice(passthrough,1,None))
            if len(new_orf) > 0:
                orfs.append(new_orf)
            else:
                exhausted = True
                
    return([''.join(orf) for orf in orfs])

def translate_rna(sequence, codon_table):
    '''Translates RNA sequence string into protein string'''
    chunk3 = lambda seq: (''.join(k) for k in zip(seq[0::3],seq[1::3],seq[2::3]))
    result = []
    for codon in chunk3(sequence): 
        result += codon_table[dna_to_rna(codon)]
    return ''.join(result)

# TODO: for testing
#assert translate_rna('ATGATGATGTCGAGCCCC') == 'MMMSHP', 'bad'

def global_alignment(a, b,match=1, mismatch=-1, gap_init=-2, gap_extend=-4):
    '''Returns a global alignment of two sequence strings using Needleman-Wunsch-Otoh algorithm.
    Default parameters are not typical; |gap_extend| should usually be less than |gap_init| 
    See Durbin et al "Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids"
    '''
    # TODO: rewrite this using generators
    # TODO: timeit, figure out why it's slow
    # TODO: read up on state machines
    # TODO: validate against BioPython, NWalign, etc.
    # TODO: representing the in-progress aligned sequences as strings is problematic because they get copied at every inset_gap(). What data structure is appropriate?
    
    inf = sys.maxsize # TODO: don't do this
    
    def maxnloc(arg):
        '''Returns (maximum, index of maximum) in a list or tuple '''
        return (max(arg), arg.index(max(arg)))
    
    def insert_gap(word, i): 
        '''inserts gap character "-" into string at position i '''
        return word[:i] + '-' + word[i:]
    
    def basematch(x,y): 
        ''' Returns score of base-base interaction'''
        # TODO: this could be replaced with similarity function like PAM or BLOSUM for protein alignment
        return match if x is y else mismatch
    
    a_aligned, b_aligned = a, b # output
    
    # scoring matrix
    M = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)]
    I = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)]
    J = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)] 
    S = [I, M, J]
    
    # traceback matrix
    BM = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)]
    BI = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)]
    BJ = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)]
    B = [BI, BM, BJ]

    # initial values
    I[0][0] = -inf
    J[0][0] = -inf
    for i in range(1, len(a)+1):
        I[i][0] = gap_init + (i-1)*gap_extend
        M[i][0] = gap_init + (i-1)*gap_extend
        J[i][0] = -inf
        J[i][1] = gap_init
    for j in range(1, len(b)+1):
        J[0][j] = gap_init + (j-1)*gap_extend
        M[0][j] = gap_init + (j-1)*gap_extend
        I[0][j] = -inf
        I[1][j] = gap_init

    # calculate scoring matrices, and put pointers in backtrack matrices
    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            # order matters when calculating matrix values
            I[i][j], BI[i][j] = maxnloc([I[i-1][j] + gap_extend,\
                                        M[i-1][j] + gap_init])
            
            J[i][j], BJ[i][j] = maxnloc([J[i][j-1] + gap_extend,\
                                        M[i][j-1] + gap_init])
            
            M[i][j], BM[i][j] = maxnloc([I[i][j],\
                                        M[i-1][j-1] + basematch(a[i-1], b[j-1]),\
                                        J[i][j]])

    # follow pointers in backtrack matrix, from bottom right to top left corner
    # TODO: rewrite this (and score calculation) to reference backtrack matrices by name instead of by B[index]
    # TODO: rewrite this using finite state machine module, if a well-supported one exists
    max_score, which_backtrack_matrix = maxnloc([I[-1][-1], M[-1][-1], J[-1][-1]])
    i,j = len(a), len(b) # our current position in the backtrack matrices
    while i>0 and j>0:
        if which_backtrack_matrix == 0:  # BI
            if BI[i][j] == 1:
                which_backtrack_matrix = 1
            i -= 1
            b_aligned = insert_gap(b_aligned, j)
        elif which_backtrack_matrix == 1:  # BM
            if BM[i][j] == 0:
                which_backtrack_matrix = 0
            elif BM[i][j] == 2:
                which_backtrack_matrix = 2
            else:
                i -= 1
                j -= 1
        elif which_backtrack_matrix == 2:  # BJ
            if BJ[i][j] == 1:
                which_backtrack_matrix = 1
            j -= 1
            a_aligned = insert_gap(a_aligned, i)
        else:
            assert False, "backtracking should never enter this state"

    # equalize the lengths
    for _ in range(i):
        b_aligned = insert_gap(b_aligned, 0)
    for _ in range(j):
        a_aligned = insert_gap(a_aligned, 0)

    return (a_aligned, b_aligned), max_score


def main():
    ''' Wraps ORF finder, translator, and sequence aligner in commandline interface '''
    # TODO: Add option to prettyprint output
    # TODO: Add option to print output in form which is convenient to pipe to other programs
    # TODO: Check input for potential mistakes
    
    # Custom codon table
    # TODO: read this from a file
    bases = ['U', 'C', 'A', 'G']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    residues =  'VVVV'+'SSSS'+'YY**'+'CC*W'+\
                'LLLL'+'PPPP'+'SSRR'+'RRRR'+\
                'IIIM'+'TTTT'+'NNKK'+'HHQQ'+\
                'FFLL'+'AAAA'+'DDEE'+'GGGG'
    codon_table = dict(zip(codons, residues))

    # TODO: for testing
    #assert codon_table['AUG'] == 'M', 'bad'
    
    parser = argparse.ArgumentParser(description='Sequence related stuff')
    parser.add_argument('file', type=argparse.FileType('r'), nargs='+',help='The sequence files. Alignment operates on the first 2 files, while everything else operates on all files.')
    parser.add_argument('--rna', help='translate DNA to RNA', dest='do_rna', default=False, action='store_true')
    parser.add_argument('--genes', help='Find ORFs and translate DNA to protein in each ORF', dest='do_genes', default=False, action='store_true')
    parser.add_argument('--align', help='globally align the first 2 sequence files given', dest='do_align', default=False, action='store_true')
    parser.add_argument('--antisense', help='operate on antisense strand instead of sense strand', dest='is_antisense', default=False, action='store_true')
    
    args = parser.parse_args()
    sequences = []
    filenames = []
    
    for f in args.file:
        if args.is_antisense:
            sequences.append(reverse_complement(f.read().upper()))
        else:
            sequences.append(f.read().upper())
        filenames.append(f.name)
        f.close()
        
    if args.do_rna:
        for f,seq in zip(filenames,sequences):
            print('From '+f+':\n',sep='\n')
            print(dna_to_rna(seq)+'\n',sep='\n')
            
    if args.do_genes:
        for f,seq in zip(filenames,sequences):
            print('From '+f+':\n',sep='\n')
            list(map(lambda w: print(*w, sep='\n'),[('DNA: '+k, 'RNA: '+dna_to_rna(k) ,'Protein: '+translate_rna(dna_to_rna(k),codon_table),'\n') for k in print_orfs(seq,codon_table)]))
            
    if args.do_align:
        if len(sequences)<2:
            print('alignment needs 2 sequences')
        else:
            (a,b), score = global_alignment(sequences[0],sequences[1])
            alignment_guide = ''.join(('|' if i is j else " " for (i,j) in zip(a,b)))
            printme = zip(a,alignment_guide,b)
            #horiz2vert = str.maketrans('-|','|-')
            #for line in printme:
            #    print(''.join(line).translate(horiz2vert))
            print(filenames[0]+':\n\t'+a+'\n\n'+filenames[1]+':\n\t'+b)
            print('\nscore: ',score)
        

if __name__ == '__main__': 
    main()

