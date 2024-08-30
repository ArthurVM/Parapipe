import os
import sys
import argparse
import pysam
import yaml
import json
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple

from Ardal import *


def dataframe_to_dict(df):
    ## convert dataframe to a dictionary
    result_dict = {}

    for allele_id, allele_data in df.items():
        allele_dict = allele_data.to_dict()
        result_dict[allele_id] = allele_dict

    return result_dict


def merge_dicts(dicts):
    merged_dict = {}

    for current_dict in dicts:
        for allele_id, allele_data in current_dict.items():
            if allele_id not in merged_dict:
                merged_dict[allele_id] = {}

            for sample_id, value in allele_data.items():
                if sample_id not in merged_dict[allele_id]:
                    merged_dict[allele_id][sample_id] = 0
                merged_dict[allele_id][sample_id] += value

    return merged_dict


def dict_to_dataframe(input_dict):
    df = pd.DataFrame.from_dict(input_dict, orient='index')
    df = df.transpose()
    df = df.fillna(0)
    return df


def mergeAMs(allele_matrix, database_matrix):
    
    dataframes = [allele_matrix, database_matrix]
    
    dicts = []
    ## Iterate through each DataFrame
    for df in dataframes:
        d = dataframe_to_dict(df)
        dicts.append(d)

    merged_dicts = merge_dicts(dicts)

    # print(dict_to_dataframe(merged_dicts))
    return dict_to_dataframe(merged_dicts)


def read_json(jfile):
    with open(jfile, 'r') as fin:
        return json.load(fin)
    

def countSNPs(ard_mat, sample_ids):

    counts_dict = defaultdict(dict)
    allale_mapping = ard_mat.toDict()

    for id in sample_ids:
        unique_alleles = ard_mat.unique([id])
        total_alleles = allale_mapping[id]
        counts_dict[id]["total_snps"] = len(total_alleles)
        counts_dict[id]["unique_snps"] = list(unique_alleles)

    with open("allele_stats.json", 'w') as fout:
        json.dump(counts_dict, fout, indent=4)


def run(args):

    ## dictionary of form { allele0 : { sample0 : 1, sample1 : 0, sample2 : 1, ... }, allele1 : { ... }, ... }
    allele_dict = defaultdict(dict)
    suffix = ".vcf.gz"

    samples = [os.path.basename(vcf).split(suffix)[0] for vcf in args.vcfs]

    samples_to_include = [] ## IDs of samples which exceed required coverage

    for vcf_path in args.vcfs:
        sample_id = os.path.basename(vcf_path).split(suffix)[0]

        fout = open(f"{sample_id}.snps.bed", 'w')
        fout.write(f"chrom start end ref alt allele_id\n")

        mapstats_json = f"{sample_id}_mapstats.json"
        jdata = read_json(mapstats_json)

        ## check coverage is sufficient
        if jdata["mapping_stats"]["coverage_breadth_hist"]["5"] >= args.mincov*100:
            samples_to_include.append(sample_id)

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

    if args.database != "false":
        database_matrix = pd.read_csv(args.database, index_col=0) 
        print(database_matrix)
        allele_matrix = mergeAMs(allele_matrix, database_matrix)
    
    allele_matrix.to_csv(f"full_allele_matrix.csv")

    ## find unique and total SNPs for all samples
    ard_mat = Ardal("./full_allele_matrix.csv")
    countSNPs(ard_mat, ard_mat.matrix.index)

    ## subset the allele matrix for only samples which exceed coverage thresholds and output
    subset_ard_mat = ard_mat.subsetbyGUID(samples_to_include)
    subset_ard_mat.matrix.to_csv("allele_matrix.csv")


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--vcfs', required=True, action='store', nargs="*", help='VCF files to contruct a merged allele matrix from')
    parser.add_argument('--qual', action='store', default=30, help='Quality score threshold for introducing a variant into a sequence.')
    parser.add_argument('--mapstats', action='store', nargs="*", help='Mapstats JSONs for quality control.')
    parser.add_argument('--database', action='store', help='Allele matrix database to merge with.')
    parser.add_argument('--mincov', action='store', default=0.8, help='The minimum fraction of the genome which must be covered to a depth of 5x to include a sample in phylogenetic analysis. Default=0.8.')

    args = parser.parse_args(argv)

    return args


if __name__=="__main__":
    args = parse_args(sys.argv)
    run(args)
