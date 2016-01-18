#!/usr/bin/python3

from Bio import SeqIO
import sys

with sys.stdin as file_handle:
    sequences = SeqIO.parse(file_handle,"fasta")
    print(sum([k.seq == k.reverse_complement().seq for k in sequences]))