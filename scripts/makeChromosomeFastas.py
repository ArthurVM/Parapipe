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

def writeFastas(chrdict, ref_id, pairwise=True):
    """ Takes a chromosome dict structured as

    {
    chr0 : { A0 : seq, A1 : seq, ... , An : seq }
    chr1 : { A0 : seq, A1 : seq, ... , An : seq }
    ... ,
    chrn : { A0 : seq, A1 : seq, ... , An : seq }
    }

    and outputs a set of multifasta files for each chromosome or a pairwise
    set to align each chromosome against its analog in the reference
    """

    for chr, subdict in chrdict.items():

        if pairwise:
            for genome_id, seq in subdict.items():
                with open(f"{genome_id}_{chr}.fasta", "w") as fout:
                    fout.write(f">{ref_id}.{chr}\n{subdict[ref_id]}\n")
                    fout.write(f">{genome_id}.{chr}\n{seq}\n")

        else:
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

    ref_path = genome_list[-1]
    ref_id = os.path.splitext(os.path.basename(ref_path))[0]
    print(f"Assuming reference is {ref_id}...")

    chrdict = defaultdict(dict)

    for genome in genome_list:
        id = os.path.splitext(os.path.basename(genome))[0]
        for chr, seq in getChrs(genome):
            chrdict[chr][id] = seq

    writeFastas(chrdict, ref_id)

if __name__=="__main__":
    if len(sys.argv) < 2:
        print("Please supply more than one assembly. USAGE: makeChromosomeFastas.py <A0> <A1> ... <An>")
        sys.exit(1)
    else:
        main(sys.argv)
