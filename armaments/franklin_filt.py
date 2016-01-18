#!/usr/bin/python3

import sys
from Bio import SeqIO

with sys.stdin as input_handle:
    # rewrite this to read from stdin
    threshold,percent = map(int,input_handle.readline().strip().split()) # threshold in first line
    sequences = SeqIO.parse(input_handle,"fastq") # fastq data in subsequent lines
    #hits = map(lambda x: mean(x)<threshold,[k.letter_annotations['phred_quality'] for k in sequences])
    #print(sum(hits))
    phred_scores = (k.letter_annotations['phred_quality'] for k in sequences)
    phred_fractions = (sum(map(lambda x:x>threshold,k))/len(k) for k in phred_scores)
    number_filtered_out = sum(map(lambda x:x<percent/100,phred_fractions))
    print(number_filtered_out+1)