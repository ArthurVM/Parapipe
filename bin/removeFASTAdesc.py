"""
Prepares a temporary FASTA file, stripping out all the sequence descriptions
"""

import sys
import os
from Bio import SeqIO

def processFasta(fasta):
    """
    Split multiFASTA into individual sequences and strips out the description.

    e.g.
    >Chr1 | chrmosome 1
    ATCGATCG...

    becomes

    >Chr1
    ATCGATCG...
    """

    fout = open(f"tmp.descStripped.fasta", "w")

    for rec in SeqIO.parse(fasta, "fasta"):
        print(f">{rec.id}\n{rec.seq}", file=fout)

    fout.close()

def main(args):
    stripdesc(args[0])

if __name__=="__main__":
    args = sys.argv[1:]
    if len(args) != 1:
        print("USAGE: removeFASTAdesc.py <FASTA>\n\nPlease provide a FASTA file!\n")
        sys.exit(1)
    else:
        main(args)
