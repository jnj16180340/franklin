#!/usr/bin/python3

from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Alphabet import generic_dna
import itertools
import sys

# Reused my ORF finder because BioPython's CDS=True requires that the sequence starts with a start codon and ends with a stop codon
# which is really irritating.
# with a little more work this could be integrated with BioPython more seamlessly
def print_orfs(seq,codon_table_ncbi):
    ''' Finds open reading frames in DNA string, including nested ORFs, and prints them, along with RNA and protein translations'''
    #stop_codons = [rna_to_dna(k) for k, v in codon_table.items() if v == '*']
    stop_codons = codon_table_ncbi.stop_codons
    #start_codons = [rna_to_dna(k) for k, v in codon_table.items() if v == 'M']
    start_codons = codon_table_ncbi.start_codons
    codon_table = codon_table_ncbi.forward_table
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

with sys.stdin as file_handle:
    temp = Seq(file_handle.read(),generic_dna)   

x = CodonTable.unambiguous_dna_by_name["Standard"]
x.start_codons = ['ATG'] # "standard" table includes 2 other alternative start codons!

nuc_orfs = sorted(print_orfs(str(temp.reverse_complement()),x)+print_orfs(str(temp),x),key=len,reverse=True)

print(str(Seq(nuc_orfs[0]).translate()))