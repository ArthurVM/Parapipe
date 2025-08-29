import os
import sys
import argparse
import yaml
import json
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple
from pysam import VariantFile


def read_json(jfile):
    with open(jfile, 'r') as fin:
        return json.load(fin)


def run(args):
    vcf_in = VariantFile(args.vcf)  # auto-detect input format
    vcf_out = VariantFile(args.out, 'w', header=vcf_in.header)

    for rec in vcf_in:

        #filter by quality
        if len(rec.alleles) > 1 and rec.qual >= args.qual:

            ## filter out indels
            if args.only_snp:
                if len(rec.alleles[0]) and len(rec.alleles[1]) == 1:
                    vcf_out.write(rec)
            else:
                vcf_out.write(rec)


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('vcf', action='store', help='VCF file to filter.')
    parser.add_argument('--qual', action='store', type=int, default=30, help='Quality score threshold.')
    parser.add_argument('--only_snp', action='store_true', default=False, help='Filter out all except SNPs.')
    parser.add_argument('-o', '--out', action='store', default="-", help='Output path.')

    args = parser.parse_args(argv)

    return args


if __name__=="__main__":
    args = parse_args(sys.argv)
    run(args)
