#!/usr/bin/python3

import sys
from Bio import SeqIO

with sys.stdin as input_handle:
    # SeqIO.parse() is an iterator. 
    # SO it keeps the file open until it's been iterated through.
    sequences = SeqIO.parse(input_handle,"fastq")
    with sys.stdout as output_handle:
        SeqIO.write(sequences,output_handle,'fasta')

