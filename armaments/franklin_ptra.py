#!/usr/bin/python3
import sys
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.Seq import translate
from operator import itemgetter

with sys.stdin as input_handle:
    dna = input_handle.readline().strip()
    protein = input_handle.readline().strip()
    is_match = [(k,translate(dna, stop_symbol='', to_stop=False, table = k) == protein) for k in CodonTable.unambiguous_dna_by_id.keys()]
    print(filter(itemgetter(1),is_match).__next__()[0])