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

def seqQC(rec, gene, assemblyID):
    """ quality control sequences, involving:
    - screening out sequences with N space
    - screening out sequences which do not start with a start codon
    - screening out sequences which to not terminate with a stop codon
    """
    nspace = rec.seq.count("N")
    if nspace >= 1:
        print(f"nspace exceeds threshold in {gene} {assemblyID} : {nspace}")
        return False

    aa_seq = rec.seq.translate()
    if aa_seq[0] != "M":
        print(f"Methionine start codon not detected in {gene} {assemblyID}")
        return False

    if aa_seq[-1] != "*":
        print(f"Stop codon not detected in {gene} {assemblyID}")
        return False

    return True

def createOrthogroups(cdslibs):

    orthodict = defaultdict(dict)
    for lib in cdslibs:
        libprefix = path.basename(lib).split("-cdna")[0]
        for rec in SeqIO.parse(lib, "fasta"):
            orthodict[rec.id][libprefix]=rec

    return orthodict

    # for lib, subdict in orthodict.items():
    #     print(lib, subdict)

def makeOrthoMultifastas(gene, libdict):

    orthofastaN = f"tmp.{gene}.fa"
    orthofastaA = f"tmp.{gene}.faa"

    ## lists for storing fasta sequences for writing to file, allowing a check of the number of sequences
    writeboxN = []
    writeboxAA = []

    with open(orthofastaN, "w") as foutN, open(orthofastaA, "w") as foutA:
        for assemblyID, rec in libdict.items():

            if not seqQC(rec, gene, assemblyID):
                continue

            ## strip out stop codon from nucleotide fasta
            recline = f">{assemblyID} | {rec.description}\n{rec.seq[:-3]}"
            writeboxN.append(recline)
            recline = f">{assemblyID} | {rec.description}\n{rec.seq.translate()}\n"
            writeboxAA.append(recline)

        if len(writeboxN)>1:
            foutN.write("\n".join(writeboxN))
            foutA.write("\n".join(writeboxAA))
            return orthofastaN

        else:
            print(f"Not enough sequences to align in {gene}. Skipping...")
            return None

def alignOrthos(orthodict):
    """ Align orthologs for dNdS analysis
    """

    outfile_box = []

    for gene, libdict in orthodict.items():
        orthofasta = makeOrthoMultifastas(gene, libdict)

        ## check there are enough genes to align
        if orthofasta == None:
            continue

        outfile = path.splitext(orthofasta)[0   ]
        outfile_box.append(outfile)

        runline = f"clustalw2 -INFILE={orthofasta} -ALIGN -TYPE=DNA -OUTPUT=FASTA -OUTPUTTREE=nexus -OUTFILE={outfile}.aln.fa"
        command(runline).run_comm(0)

    return outfile_box

def gen_argparser(argv):

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)

    parser.add_argument('gene_libs', type=is_file, action='store', nargs="*",
     help='gene sequences in FASTA format.')

    args = parser.parse_args(argv)

    return args

def main(argv):

    args = gen_argparser(argv)

    orthodict = createOrthogroups(args.gene_libs)
    outfile_box = alignOrthos(orthodict)

if __name__=="__main__":
    main(sys.argv)
