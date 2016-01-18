#!/usr/bin/python3
import sys

from Bio import Entrez
from Bio import SeqIO
#from functools import reduce

def main(argv):
    Entrez.email = "jnj.16180340@rosalind.info"
    # input() reads stdin
    #query_string = reduce(lambda x,y: x+', '+y, input().strip().split())
    query_string = input().strip().split()
    
    handle = Entrez.efetch(db="nucleotide", id=query_string, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta")) #we get the list of SeqIO objects in FASTA format
    records.sort(key=lambda x:len(x.seq),reverse=True)
    print(records.pop().format('fasta'))

if __name__ == "__main__":
    main(sys.argv)