import os
import sys
import argparse
import pysam
import yaml
import json
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple


def read_json(jfile):
    with open(jfile, 'r') as fin:
        return json.load(fin)


def run(args):

    ## dictionary of form { allele0 : { sample0 : 1, sample1 : 0, sample2 : 1, ... }, allele1 : { ... }, ... }
    allele_dict = defaultdict(dict)
    suffix = ".vcf.gz"

    samples = [os.path.basename(vcf).split(suffix)[0] for vcf in args.vcfs]

    for vcf_path in args.vcfs:
        sample_id = os.path.basename(vcf_path).split(suffix)[0]

        fout = open(f"{sample_id}.snps.bed", 'w')
        fout.write(f"chrom start end ref alt allele_id\n")

        mapstats_json = f"{sample_id}_mapstats.json"
        jdata = read_json(mapstats_json)

        ## check at least 95% of the genome is covered to at least 5x
        ## if not, exclude this sample from SNP phylo tree construction
        # if jdata["mapping_stats"]["coverage_breadth_hist"]["5"] < 90.0:
        #     print(f"ignored {vcf_path}")
        #     continue

        vcf = pysam.VariantFile(vcf_path)

        for record in vcf.fetch():
            ## filter for high quality SNPs only
            if record.qual >= args.qual and len(record.alleles[1])==1:
                allele_id = f"{record.chrom}.{record.pos}.{record.alleles[1]}"

                fout.write(f"{record.chrom} {record.pos} {record.pos+1} {record.alleles[0]} {record.alleles[1]} {allele_id}\n")

                if allele_id not in allele_dict:
                    for s in samples:
                        allele_dict[allele_id][s] = 0

                allele_dict[allele_id][sample_id] = 1

        fout.close()

    allele_matrix = pd.DataFrame.from_dict(allele_dict)
    allele_matrix.to_csv(f"allele_matrix.csv")


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--vcfs', required=True, action='store', nargs="*", help='VCF files to contruct a merged allele matrix from')
    parser.add_argument('--qual', action='store', default=30, help='Quality score threshold for introducing a variant into a sequence.')
    parser.add_argument('--mapstats', action='store', nargs="*", help='Mapstats JSONs for quality control.')

    args = parser.parse_args(argv)

    return args


if __name__=="__main__":
    args = parse_args(sys.argv)
    run(args)
