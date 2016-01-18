#!/usr/bin/python3

import sys
from statistics import mean
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools

output_sequences = []

with sys.stdin as input_handle:
    threshold = int(input_handle.readline().strip()) # threshold in first line
    sequences = SeqIO.parse(input_handle,"fastq") # fastq data in subsequent lines
    for k in sequences:
        # TODO: this is a mess
        base_quality = zip(k.seq,k.letter_annotations['phred_quality'])
        base_quality = itertools.dropwhile(lambda l: l[1]<threshold,base_quality)
        base_quality = reversed(list(base_quality))
        base_quality = itertools.dropwhile(lambda l: l[1]<threshold,base_quality)
        base_quality = reversed(list(base_quality))
        base_quality = list(base_quality)
        
        new_seq = ''.join([p[0] for p in base_quality])
        new_lett_anno = {'phred_quality':[p[1] for p in base_quality]}
        
        # copy each sequence record, but replace seq and letter_annotations with trimmed versions
        output_sequences.append(SeqRecord(new_seq, id=k.id, name=k.name, description=k.description, dbxrefs=k.dbxrefs, features=k.features, annotations=k.annotations, letter_annotations=new_lett_anno))
        
SeqIO.write(output_sequences,sys.stdout,'fastq')