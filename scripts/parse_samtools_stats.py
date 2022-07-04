#!/usr/bin/env python3

from os import path
import argparse
import json
import pysam
import numpy as np


def total_genome_size_from_bam(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    try:
        total_length = sum(bam.lengths)
    except:
        raise Exception(
            f"Error getting total of reference seq lengths from BAM file {bam_file}"
        )

    bam.close()
    return total_length


def coverage_hist_to_coverage_breadth(coverage_hist, genome_size, max_cov=100):
    cumulative_cov = {max_cov: 0}
    total_bases = 0
    for cov, depth in sorted(coverage_hist.items(), reverse=True):
        if cov >= max_cov:
            cumulative_cov[max_cov] += depth
        else:
            cumulative_cov[cov] = depth

        total_bases += depth

    depths = [1, 2, 3, 4, 5, 10, 15, 25, 50, 100]
    cov_breadth = {
        d: sum([v for k, v in cumulative_cov.items() if k >= d]) for d in depths
    }
    return {k: round(100 * v / genome_size, 2) for k, v in cov_breadth.items()}


def parse_samtools_stats_file(infile, doc_bed, genome_size):
    stats = {}
    wanted_keys = {
        "raw total sequences:": False,
        "reads mapped:": False,
        "reads duplicated:": False,
        "bases mapped (cigar):": False,
        "bases trimmed:": False,
        "error rate:": True,
        "average quality:": True,
        "insert size average:": True,
        "insert size standard deviation:": True,
        "inward oriented pairs:": False,
        "outward oriented pairs:": False,
        "pairs with other orientation:": False,
    }

    coverage_hist = {}

    with open(infile) as f:
        for line in f:
            if line.startswith("SN"):
                fields = line.rstrip().split("\t")
                if fields[1] in wanted_keys:
                    key = (
                        fields[1]
                        .replace(" ", "_")
                        .rstrip(":")
                        .replace("(", "")
                        .replace(")", "")
                    )
                    value = (
                        float(fields[2]) if wanted_keys[fields[1]] else int(fields[2])
                    )

                    stats[key] = value
            elif line.startswith("COV"):
                try:
                    _, _, pos, count = line.rstrip().split("\t")
                    pos = int(pos)
                    count = int(count)
                except:
                    raise Exception(
                        f"Error parsing this COV line in file {infile}: {line}"
                    )
                coverage_hist[pos] = count

    mean_DOC, median_DOC = parse_doc_bed(doc_bed)
    stats["mean_depth_of_coverage"] = mean_DOC
    stats["median_depth_of_coverage"] = median_DOC

    stats["coverage_breadth"] = coverage_hist_to_coverage_breadth(
        coverage_hist, genome_size
    )
    stats["genome_size"] = genome_size
    return stats

def calc_GG_area(args, stats):

    ## Calculates the area under the GG curve

    with open(args.gg_file, "r") as fin:
        ggdict = {int(line.split("\t")[0]) : float(line.split("\t")[1]) for line in fin.readlines()}

    gg_array = np.array([float(g) for w, g in ggdict.items()])
    window_array = np.array([int(w) for w, g in ggdict.items()])
    w = window_array[1]-window_array[0]
    norm_array = [x+1-np.max(gg_array) for x in gg_array]

    m_area = np.trapz(gg_array, dx=w)/np.max(window_array)
    m_areaN = np.trapz(norm_array, dx=w)/np.max(window_array)
    # print(f"GG-a={np.round(m_area, 3)}\tnGG-a={np.round(m_areaN, 3)}")

    stats["gg_array"] = ggdict
    stats["GG_area"] = np.round(m_area, 3)
    stats["nGG_area"] = np.round(m_areaN, 3)
    return stats

def parse_doc_bed(doc_bed):

    with open(doc_bed, 'r') as fin:
        cov_list = [int(line.strip('\n').split('\t')[-1]) for line in fin.readlines()]

    mean_DOC = np.mean(cov_list)
    median_DOC = np.median(cov_list)

    return mean_DOC, median_DOC

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Gather QC stats from BAM file and output from samtools stats, prints JSON to stdout",
        usage="%(prog)s <bam> <samtools_stats> <GG>",
    )
    parser.add_argument("bam", help="sorted indexed BAM file")
    parser.add_argument("samtools_stats", help="File made by samtools stats bam")
    parser.add_argument("doc_bed", help="Depth of coverage BED file generated using Samtools depth")
    parser.add_argument("gg_file", help="File containing gini stats")
    options = parser.parse_args()
    genome_size = total_genome_size_from_bam(options.bam)
    stats = parse_samtools_stats_file(options.samtools_stats, options.doc_bed, genome_size)
    stats = calc_GG_area(options, stats)
    tmp = dict()
    tmp.update({"bam_file" : path.basename(options.bam)})
    tmp.update({"mapping_stats":stats})
    print(json.dumps(tmp, sort_keys=True, indent=2))
