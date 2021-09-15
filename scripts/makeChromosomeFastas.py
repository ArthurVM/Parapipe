#!/usr/bin/env python3

import sys
import os
import re
from collections import defaultdict
from Bio import SeqIO

def getChrs(fasta_file):
    """ takes a fasta file and yields chromosomes
    """
    for rec in SeqIO.parse(fasta_file, "fasta"):
        yield rec.id, rec.seq

def writeFastas(chrdict):
    """ Takes a chromosome dict structured as

    {
    chr0 : { A0 : seq, A1 : seq, ... , An : seq }
    chr1 : { A0 : seq, A1 : seq, ... , An : seq }
    ... ,
    chrn : { A0 : seq, A1 : seq, ... , An : seq }
    }

    and outputs a set of multifasta files for each chromosome.
    """
    for chr, subdict in chrdict.items():
        with open(f"{chr}.fasta", "w") as fout:
            for genome_id, seq in subdict.items():
                fout.write(f">{genome_id}.{chr}\n{seq}\n")

def parseArgs(args):
    """ reads arguments and attempts to interpret as a literal, returns a list of assemblies
    """
    # TODO: God this is so horrible, need to find a nicer way to parse nextflow lists as arguments
    if args[1].startswith("["):
        return [i for i in re.split("[][, ]", args[1]) if i != ""]
    else:
        return args[1:]

def main(args):
    """ main function
    """
    genome_list = parseArgs(args)

    chrdict = defaultdict(dict)

    for genome in genome_list:
        id = os.path.splitext(os.path.basename(genome))[0]
        for chr, seq in getChrs(genome):
            chrdict[chr][id] = seq

    writeFastas(chrdict)

if __name__=="__main__":
    if len(sys.argv) < 2:
        print("Please supply more than one assembly. USAGE: makeChromosomeFastas.py <A0> <A1> ... <An>")
        sys.exit(1)
    else:
        main(sys.argv)
