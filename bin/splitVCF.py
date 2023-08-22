"""
Prepares temporary vcf files, splitting by chromosome
"""

import sys
import os
from collections import defaultdict

def splitVCF(vcf):
    """ Reads a vcf file and splits it into seperate files according to chromosome ID
    """

    meta = []
    chrdict = defaultdict(list)
    with open(vcf, "r") as fin:
        for line in fin.readlines():
            line = line.strip("\n")
            if line.startswith("#"):
                meta.append(line)
            else:
                chr = line.split("\t")[0]
                chrdict[chr].append(line)

    outdir = f"{os.getcwd()}/vcf_split"
    os.mkdir(outdir)

    for chr, line in chrdict.items():
        with open(f"{outdir}/tmp.{chr}.vcf", "w") as fout:
            newmeta = []
            for m in meta:
                if m.startswith("##contig"):
                    if chr in m:
                        newmeta.append(m)
                    else:
                        continue
                else:
                    newmeta.append(m)

            fout.write("\n".join(newmeta) + "\n")
            fout.write("\n".join(line))

def main(args):
    splitVCF(args[0])

if __name__=="__main__":
    args = sys.argv[1:]
    if len(args) != 1:
        print("USAGE: splitVCF.py <VCF>\n\nPlease provide a VCF file!\n")
        sys.exit(1)
    else:
        main(args)
