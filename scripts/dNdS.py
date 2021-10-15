#!/usr/bin/env python3

"""
aligns orthogroups
"""

# TODO: EMBOSS water reader class

import sys
import argparse
import subprocess
from utilities import *
from collections import defaultdict
from os import path, listdir, walk
from Bio import SeqIO
from Bio.Seq import Seq

def readalnfile(aln_file):
    return [ rec for rec in SeqIO.parse(aln_file, "fasta") ]

def processSeq(seq):
    ## performs sequence QC and processes for comparison
    # offcut = len(seq)%3
    # if offcut != 0:
    #     print(f"WARNING: Sequence length not a multiple of 3, trimming by {offcut}...")
    #     seq = seq[:-offcut]
    #
    # codons = chunk(seq)
    # aa_seq = [str(c.translate()) if "-" not in c else str(c) for c in codons]
    # print(aa_seq, map(lambda x: str(x), codons))

    if seq[0] != "M":
        print(f"WARNING: Methionine start codon not detected.")

def readFastaFile(aln_file):
    ## takes an aln_file and reads in the FASTA file containing nucleotide sequence

    prefix = aln_file.split(".aln.fa")[0]
    return [rec for rec in SeqIO.parse(f"{prefix}.fa", "fasta")]

def getCodons(aa_seq, n_seq):
    ## creates a list of tuples containing the amino acid and nucleotide sequence for an amino acid sequence
    seq_tuple_box = []

    c = 0
    for aa in aa_seq:

        ## skip gaps in the alignment
        if aa == "-":
            seq_tuple_box.append(["-", "-"])
            continue

        codon = n_seq[3*c:3*c+3]
        # sys.stdout.write(f"{c} {codon}\t")

        ## check that nothing has gone wrong with pairing
        # sys.stdout.write(f"{str(Seq(codon).translate())} == {aa}\t")

        if not str(Seq(codon).translate()) == aa:
            sys.stdout.write(f"{str(Seq(codon).translate())} == {aa}\t")

        seq_tuple_box.append([aa, codon])
        c+=1

    return seq_tuple_box


def checkSites(aa_pairs, ref_seq, q_seq):
    """Get the sum of synonymous sites from a zipped amino acid alignment"""
    syn_site = 0
    non_syn_subs = 0
    ambiguous_sites = 0
    for pair in aa_pairs:
        if "X" in pair:
            ambiguous_sites+=1
            continue
        elif set(pair)==1:
            non_syn_subs+=1
        else:
            syn_site+=1



def dNdSmain(aln_file):
    print(aln_file)

    fasta_n_box = readFastaFile(aln_file)
    alnbox = readalnfile(aln_file)

    processSeq(alnbox[0].seq)

    seq_tuple_dict = { alnbox[i].id : getCodons(str(alnbox[i].seq), str(fasta_n_box[i].seq)) for i in range(len(alnbox)) }

    # print(seq_tuple_dict)

    # ref=alnbox[0]
    # for rec in alnbox[1:]:
    #     print(f"\t{rec.seq}")
    #     aa_pairs = zip(ref.seq, rec.seq)
    #     ambiguous_sites(aa_pairs, )

def gen_argparser(argv):

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)

    parser.add_argument('aln_dir', type=is_dir, action='store',
     help='directory containing alignment files.')

    args = parser.parse_args(argv)

    return args

def main(argv):

    args = gen_argparser(argv)

    alnfile_box = [f for f in listdir(args.aln_dir) if f.endswith(".aln.fa")]

    for aln_file in alnfile_box:
        dNdSmain(aln_file)

if __name__=="__main__":
    main(sys.argv)
