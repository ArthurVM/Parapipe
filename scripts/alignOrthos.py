#!/usr/bin/env python3

"""
aligns orthogroups
"""

# TODO: EMBOSS water reader class

import sys
import argparse
import subprocess
from collections import defaultdict
from os import path, listdir, walk
from Bio import SeqIO

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

    orthofasta = f"tmp.{gene}.fa"
    with open(orthofasta, "w") as fout:
        for assemblyID, rec in libdict.items():
            recline = f">{assemblyID} | {rec.description}\n{rec.seq}\n"
            fout.write(recline)

    return orthofasta

def alignOrthos(orthodict):
    """ Align orthologs for dNdS analysis
    """

    outfile_box = []

    for gene, libdict in orthodict.items():
        orthofasta = makeOrthoMultifastas(gene, libdict)

        outfile = path.splitext(orthofasta)[0]
        outfile_box.append(outfile)

        runline = f"clustalw2 -INFILE={orthofasta} -ALIGN -TYPE=DNA -OUTPUT=FASTA -OUTFILE={outfile}.aln.fa"
        subprocess.run(runline, shell=True, check=True)

    return outfile_box

def is_file(filename):
    """ Checks if a path is a file """

    if not path.isfile(filename):
        time=strftime("%H:%M:%S", localtime())
        print("{time} :: No file found at {f}".format(
                time=time,
                f=filename), end="\n", file=sys.stderr, flush=True)
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))

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
