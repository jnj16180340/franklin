#!/usr/bin/python3
import string
import sys

from Bio.Seq import Seq

def unpunctuate(s):
    return s.translate(str.maketrans("","", string.punctuation))

def main(argv):
    # input() reads stdin
    dataseq = Seq(input().strip())
    print(unpunctuate(str([dataseq.count(k) for k in 'ACGT'])))

if __name__ == "__main__":
    main(sys.argv)