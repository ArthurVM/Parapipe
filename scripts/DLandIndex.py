"""
Downloads genome sequence and annotation files and prepares index files for mapping.
"""

import sys
import os
import argparse
from utilities import *

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message

        sys.stderr.write(f"INPUT ERROR : {self.expression}\n{self.message}")
        sys.exit(1)

def checkGenomeID(args):
    """ checks a genome id exists within the list of accepted genomes, throws and error and exits upon failure
    """
    genome_id = args.genome_id
    scriptdir = os.path.split(args.script_path)[0]
    with open(f"{scriptdir}/../resources/genome_list.txt", "r") as fin:
        accepted_ids = [line.strip("\n") for line in fin.readlines()]

    if genome_id not in accepted_ids:
        ## throw an error if the genome doesn't exist in genome_list.txt and suggest similar genomes

        tmp_id = genome_id.lower()
        similar_ids = [id for id in accepted_ids if id.split("_")[0] == tmp_id.split("_")[0]]

        if tmp_id in accepted_ids:
            InputError(f"checkGenomeID", f"'{genome_id}' is not an accepted genome id. Did you mean {tmp_id}?\n")
        elif len(similar_ids) != 0:
            jids = ", ".join(similar_ids)
            InputError(f"checkGenomeID", f"'{genome_id}' is not an accepted genome id. Did you mean any of the following: {jids}?\n")
        else:
            InputError(f"checkGenomeID", f"'{genome_id}' is not an accepted genome id.\n")

    else:
        return genome_id.lower()

def downloadData(genome_id):
    """ Downloads genome and annotation data for the specified reference genome
    """

    if genome_id.startswith("cryptosporidium"):
        ## handle all Crypto seperately

        ftpbase = "https://cryptodb.org/common/downloads/Current_Release/"
        cryptorefdict = {\
        "hominis" : "ChominisUdeA01", \
        "parvum" : "CparvumIOWA-ATCC", \
        "muris" : "CmurisRN66", \
        "meleagridis" : "CmeleagridisUKMEL1"
        }
        species = genome_id.split("_")[1]
        ftpdir = os.path.join(ftpbase, cryptorefdict[species])

        gff_path = getCryptoGFF(ftpdir)
        gaf_path = getCryptoGAF(ftpdir)
        fasta_path, annoCDS_path = getCryptoFASTA(ftpdir)

        print(gff_path, gaf_path, fasta_path, annoCDS_path)

        ## perform downloads
        command(f"curl {gff_path} -o ./{genome_id}.gff").run_comm(0)
        command(f"curl {gaf_path} -o ./{genome_id}_GO.gaf").run_comm(0)
        command(f"curl {fasta_path} -o ./{genome_id}.fasta").run_comm(0)
        command(f"curl {annoCDS_path} -o ./{genome_id}_cds.fasta").run_comm(0)

    else:
        ftpdir = "ftp://ftp.ensemblgenomes.org/pub/current/protists/fasta/"
        None

def getCryptoGFF(ftpdir):
    """ takes html output from the curl runline of CryptoDB and returns the GFF path
    """
    cmdline = f"curl -s {ftpdir}/gff/data/"
    curlOut = str(command(cmdline).run_comm(1).decode("utf-8").rstrip())
    gff = curlOut.strip("\n").split('.gff">')[1].split('</a>')[0]

    return os.path.join(f"{ftpdir}/gff/data/", gff)

def getCryptoGAF(ftpdir):
    """ takes html output from the curl runline of CryptoDB and returns the GAF path
    """
    cmdline = f"curl -s {ftpdir}/gaf/"
    curlOut = str(command(cmdline).run_comm(1).decode("utf-8").rstrip())
    gff = curlOut.strip("\n").split('.gaf">')[1].split('</a>')[0]

    return os.path.join(f"{ftpdir}/gaf/data/", gff)

def getCryptoFASTA(ftpdir):
    """ takes html output from the curl runline of CryptoDB and returns the FASTA path
    """
    cmdline = f"curl -s {ftpdir}/fasta/data/"
    curlOut = str(command(cmdline).run_comm(1).decode("utf-8").rstrip())
    fasta =  curlOut.strip("\n").split('Genome.fasta">')[1].split('</a>')[0]
    annotatedCDS = curlOut.strip("\n").split('AnnotatedCDSs.fasta">')[1].split('</a>')[0]

    return os.path.join(f"{ftpdir}/fasta/data/", fasta), os.path.join(f"{ftpdir}/fasta/data/", annotatedCDS)

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('genome_id', action='store', help='A reference genome ID from the genome_list.txt file.')

    args = parser.parse_args(argv)
    return args

def main(argv):

    args = parseArgs(argv)
    genome_id = checkGenomeID(args)
    downloadData(genome_id)

if __name__=="__main__":
    main(sys.argv)
