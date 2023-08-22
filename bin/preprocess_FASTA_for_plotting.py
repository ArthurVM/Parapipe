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

    for rec in SeqIO.parse(fasta, "fasta"):
        fout = open(f"tmp.{rec.id}.fasta", "w")
        print(f">{rec.id}\n{rec.seq}", file=fout)
        fout.close()

def main(args):
    processFasta(args[0])

if __name__=="__main__":
    args = sys.argv[1:]
    if len(args) != 1:
        print("USAGE: preprocess_FASTA_for_plotting.py <FASTA>\n\nPlease provide a FASTA file!\n")
        sys.exit(1)
    else:
        main(args)
