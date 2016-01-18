#!/usr/bin/python3

from Bio import Entrez
from Bio import SeqIO
import subprocess
import sys
import re

def main(argv):
    Entrez.email = 'jnj.16180340@rosalind.info'
    query_strings = input().strip().split()
    
    for query in query_strings:
        with open(query+'.gbk','w') as out_handle:
            # Entrez.efetch is not trivial to contextualize
            in_handle = Entrez.efetch(db='nucleotide',id=query,rettype='gb',retmode='text')
            out_handle.write(in_handle.read())
            in_handle.close()
       
    stretcher_command = 'stretcher -auto -snucleotide1 -snucleotide2 -datafile EDNAFULL -gapopen 10 -gapextend 1 -aformat3 score -aname3 rosalind_need -aextension3 stretcher -asequence '+query_strings[0]+'.gbk -bsequence '+query_strings[1]+'.gbk '
    subprocess.check_call(stretcher_command,shell=True)
    
    with open('rosalind_need.stretcher','r') as aln_handle:
        # use http://pythex.org/
        # \((.{0,})\)
        # \((.*)\)
        rex = '\((.*)\)'
        print(re.findall(rex,aln_handle.read())[0])

if __name__ == "__main__":
    main(sys.argv)
    subprocess.check_call('rm *.gbk *.stretcher',shell=True)