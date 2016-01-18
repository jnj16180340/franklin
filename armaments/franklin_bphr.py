#!/usr/bin/python3

import sys
from statistics import mean
from Bio import SeqIO

with sys.stdin as input_handle:
    threshold = int(input_handle.readline().strip()) # threshold in first line
    sequences = SeqIO.parse(input_handle,"fastq") # fastq data in subsequent lines
    qualities = [k.letter_annotations['phred_quality'] for k in sequences]
    
print(sum(list(map(lambda q: mean(list(q))<threshold, zip(*qualities)))))