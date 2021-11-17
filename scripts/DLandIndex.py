"""
Downloads genome sequence and annotation files and prepares index files for mapping.
"""

import sys
import os
import argparse
from utilities import *
from shutil import copyfile

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

class RefData(object):
    """ Base class for parsing reference data for download """

    def __init__(self, ref_id):
        """ Reads in data from genome_list.txt and assigns variables """
        scriptdir = sys.path[0]
        with open(f"{scriptdir}/../resources/genome_list.txt") as fin:
            refdata_lines = [line.strip("\n").split("\t") for line in fin.readlines()]
            source, ref_id, path, geneID = checkGenomeID(ref_id, refdata_lines)

            self.source = source
            self.id = ref_id
            self.path = path
            self.geneID = geneID

def checkGenomeID(ref_id, refdata_lines):
    """ checks a genome id exists within the list of accepted genomes, throws and error and exits upon failure
    """
    accepted_ids = [i[1] for i in refdata_lines]

    if ref_id not in accepted_ids:
        ## throw an error if the genome doesn't exist in genome_list.txt and suggest similar genomes

        tmp_id = ref_id.lower()
        similar_ids = [id for id in accepted_ids if id.split("_")[0] == tmp_id.split("_")[0]]

        if tmp_id in accepted_ids:
            InputError(f"checkGenomeID", f"'{ref_id}' is not an accepted genome id. Did you mean {tmp_id}?\n")
        elif len(similar_ids) != 0:
            jids = ", ".join(similar_ids)
            InputError(f"checkGenomeID", f"'{ref_id}' is not an accepted genome id. Did you mean any of the following: {jids}?\n")
        else:
            InputError(f"checkGenomeID", f"'{ref_id}' is not an accepted genome id.\n")

    else:
        for refline in refdata_lines:
            if refline[1].lower() == ref_id.lower():
                return refline

def downloadData_OLD(genome_id):
    """ Downloads genome and annotation data for the specified reference genome
    """

    if genome_id.startswith("cryptosporidium"):
        ## handle all Crypto seperately

        ftpbase = "https://cryptodb.org/common/downloads/Current_Release/"
        cryptorefdict = {\
        "hominis"     : "ChominisUdeA01", \
        "parvum"      : "CparvumIOWA-ATCC", \
        "muris"       : "CmurisRN66", \
        "meleagridis" : "CmeleagridisUKMEL1"
        }
        species = genome_id.split("_")[1]
        ftpdir = os.path.join(ftpbase, cryptorefdict[species])

        gff_path = getGFF(ftpdir)
        gaf_path = getGAF(ftpdir)
        fasta_path, annoCDS_path = getFASTA(ftpdir)

        print(gff_path, gaf_path, fasta_path, annoCDS_path)

        ## perform downloads
        command(f"curl {gff_path} -o ./{genome_id}.gff").run_comm(0)
        command(f"curl {gaf_path} -o ./{genome_id}_GO.gaf").run_comm(0)
        command(f"curl {fasta_path} -o ./{genome_id}.fasta").run_comm(0)
        command(f"curl {annoCDS_path} -o ./{genome_id}_cds.fasta").run_comm(0)

    elif genome_id.startswith("plasmodium"):
        ## handle all Plasmo seperately

        ftpbase = "https://plasmodb.org/common/downloads/Current_Release/"
        plasmorefdict = {\
        "falciparum"      : "Pfalciparum3D7", \
        "falciparum_3D7"  : "Pfalciparum3D7", \
        "falciparum_7G8"  : "Pfalciparum7G8", \
        "falciparum_CD01" : "PfalciparumCD01", \
        "falciparum_Dd2"  : "PfalciparumDd2", \
        "falciparum_GA01" : "PfalciparumGA01", \
        "falciparum_GB4"  : "PfalciparumGB4", \
        "falciparum_GN01" : "PfalciparumGN01", \
        "falciparum_HB3"  : "PfalciparumHB3", \
        "falciparum_IT"   : "PfalciparumIT", \
        "falciparum_KE01" : "PfalciparumKE01", \
        "falciparum_KH01" : "PfalciparumKH01", \
        "falciparum_KH02" : "PfalciparumKH02", \
        "falciparum_ML01" : "PfalciparumML01", \
        "falciparum_SD01" : "PfalciparumSD01", \
        "falciparum_SN01" : "PfalciparumSN01", \
        "falciparum_TG01" : "PfalciparumTG01", \
        "knowlesi"        : "PknowlesiA1H1", \
        "malariae"        : "PmalariaeUG01", \
        "vivax"           : "PvivaxP01"
        }
        species = genome_id.split("_")[1]
        ftpdir = os.path.join(ftpbase, plasmorefdict[species])

        gff_path = getGFF(ftpdir)
        gaf_path = getGAF(ftpdir)
        fasta_path, annoCDS_path = getFASTA(ftpdir)

        print(gff_path, gaf_path, fasta_path, annoCDS_path)

        ## perform downloads
        command(f"curl {gff_path} -o ./{genome_id}.gff").run_comm(0)
        command(f"curl {gaf_path} -o ./{genome_id}_GO.gaf").run_comm(0)
        command(f"curl {fasta_path} -o ./{genome_id}.fasta").run_comm(0)
        command(f"curl {annoCDS_path} -o ./{genome_id}_cds.fasta").run_comm(0)

    else:
        ftpdir = "ftp://ftp.ensemblgenomes.org/pub/current/protists/fasta/"
        None

def getData(genome_id):
    """ Downloads genome and annotation data for the specified reference genome
    """

    ref = RefData(genome_id)

    if ref.source == "LOCAL":
        ## deal with local reference files
        gff_path = os.path.join(ref.path, ref.id + ".gff")
        gaf_path = os.path.join(ref.path, ref.id + "_GO.gaf")
        fasta_path = os.path.join(ref.path, ref.id + ".fasta")
        annoCDS_path = os.path.join(ref.path, ref.id + "_cds.fasta")

        print(gff_path, gaf_path, fasta_path, annoCDS_path)
        cwd = os.getcwd()
        copyfile(gff_path, os.path.join(cwd, ref.id + ".gff"))
        copyfile(gaf_path, os.path.join(cwd, ref.id + "_GO.gaf"))
        copyfile(fasta_path, os.path.join(cwd, ref.id + ".fasta"))
        copyfile(annoCDS_path, os.path.join(cwd, ref.id + "_cds.fasta"))

    elif ref.source == "FTP":
        gff_path = getGFF(ref.path)
        gaf_path = getGAF(ref.path)
        fasta_path, annoCDS_path = getFASTA(ref.path)

        print(gff_path, gaf_path, fasta_path, annoCDS_path)

        ## perform downloads
        command(f"curl {gff_path} -o ./{genome_id}.gff").run_comm(0)
        command(f"curl {gaf_path} -o ./{genome_id}_GO.gaf").run_comm(0)
        command(f"curl {fasta_path} -o ./{genome_id}.fasta").run_comm(0)
        command(f"curl {annoCDS_path} -o ./{genome_id}_cds.fasta").run_comm(0)

    else:
        InputError(f"getData", f"'{ref.source}' is not an accepted source location. Please specify FTP or LOCAL.\n")

def getGFF(ftpdir):
    """ takes html output from the curl runline of CryptoDB and returns the GFF path
    """
    cmdline = f"curl -s {ftpdir}/gff/data/"
    curlOut = str(command(cmdline).run_comm(1).decode("utf-8").rstrip())
    gff = curlOut.strip("\n").split('.gff">')[1].split('</a>')[0]

    return os.path.join(f"{ftpdir}/gff/data/", gff)

def getGAF(ftpdir):
    """ takes html output from the curl runline of CryptoDB and returns the GAF path
    """
    cmdline = f"curl -s {ftpdir}/gaf/"
    curlOut = str(command(cmdline).run_comm(1).decode("utf-8").rstrip())
    gaf = curlOut.strip("\n").split('.gaf">')[1].split('</a>')[0]

    return os.path.join(f"{ftpdir}/gaf/", gaf)

def getFASTA(ftpdir):
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
    getData(args.genome_id)

if __name__=="__main__":
    main(sys.argv)
