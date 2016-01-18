#!/usr/bin/python3
import sys
import xml.etree.ElementTree as ElementTree
import subprocess

# TODO this doesn't always work... pick better params for meme
# would prefer running this on the meme server, but let's use local copy for now
#  cat rosalind_meme.txt | ./meme_4.11.0/bin/meme -protein -mod oops -nmotifs 1 -minsites 20 -nostatus stdin

def main(argv):
    subprocess.call('rm ./meme_out/*',shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    meme_query = './meme_4.11.0/bin/meme -protein -mod oops -nmotifs 1 -minsites 20 -nostatus stdin'
    # this produces warnings because the meme xml->html converter doesn't work
    subprocess.call('cat rosalind_meme.txt | '+meme_query,shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    tree = ElementTree.parse('./meme_out/meme.xml')
    # root = ElementTree.fromstring(input().strip())
    root = tree.getroot()
    result = root.find('.//regular_expression')
    print(result.text.translate(str.maketrans('','','\n')))
    
if __name__ == "__main__":
    main(sys.argv)