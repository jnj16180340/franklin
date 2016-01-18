#!/usr/bin/python3
#import string
import sys
#import itertools
#import operator

from Bio import ExPASy
from Bio import SwissProt

def main(argv):
    # input() reads stdin
    handle = ExPASy.get_sprot_raw(input().strip()) #you can give several IDs separated by commas
    record = SwissProt.read(handle) # use SwissProt.parse for multiple proteins
    
    # there ought to be a better way to pull GO information from the record! maybe there is...
    for p in filter(lambda x:x[0]=='GO' and x[2].startswith('P:'),record.cross_references):
        print(p[2][2:])

if __name__ == "__main__":
    main(sys.argv)