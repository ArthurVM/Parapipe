import os
import shutil
import sys
import argparse
import subprocess
import pysam
import yaml
import json
import numpy as np
from collections import defaultdict, namedtuple
from Bio import SeqIO

NV = namedtuple("NV", ["chr", "start", "end", "id"])


def make_multifasta(profile_id, concats):

    out_path = f"{profile_id}.concats.fasta"
    with open(out_path, 'w') as fout:
        for fa in concats:
            for rec in SeqIO.parse(fa, "fasta"):
                fout.write(f">{rec.id}\n{rec.seq}\n")

    return out_path


def read_yaml(yaml_file):

    with open(yaml_file, 'r') as fin:
        return yaml.safe_load(fin)


def read_bed(bed):
    locdict = defaultdict(NV)

    with open(bed, 'r') as fin:
        for line in fin.readlines():
            chr, s, e, id = line.strip("\n").split(" ")
            # id = id.split("=")[1]
            locdict[id] = NV(chr, int(s), int(e), id)

    return locdict


def aln_and_tree(var_concats):
    """ Requires clustalw and iqtree2 in PATH
    """

    for profile_id, multifasta in var_concats.items():
        # out_path = make_multifasta(profile_id, concats)

        aln_line = f"clustalw2 -INFILE={multifasta} -ALIGN -OUTPUT=nexus -TYPE=DNA -OUTFILE={profile_id}.nxs"
        subprocess.call(aln_line, shell=True)

        tree_line = f"iqtree2 -s {profile_id}.nxs -b 100"
        subprocess.call(tree_line, shell=True)


def make_JSON(var_profiles):

    json_dict = defaultdict(dict)

    for profile_id, bed in var_profiles.items():
        json_dict[profile_id]["bed"] = bed

        contree = profile_id + ".nxs.contree"
        iqtree_log = profile_id + ".nxs.iqtree"
        mldist = profile_id + ".nxs.mldist"
        nxs = profile_id + ".nxs"

        with open(contree, 'r') as fin:
            json_dict[profile_id]["contree"] = fin.read().strip("\n")

        with open(mldist, 'r') as fin:
            mlbox = []
            for line in fin.readlines():
                sline = [l for l in line.strip("\n").split(" ") if l!='']
                mlbox.append(sline)

            json_dict[profile_id]["mldist"] = mlbox

    with open("iqtree.json", 'w') as fout:
        json.dump(json_dict, fout, indent=4)


def get_locus_seq(ref, scheme, id, nv, bam_path, prefix, fout):

    cov_warnings = []

    bam_index_line = f"samtools index {bam_path}"
    print(bam_index_line, flush=True)
    subprocess.call(bam_index_line, shell=True)

    cov_line = f"samtools depth {bam_path} -r \"{nv.chr}:{nv.start}-{nv.end}\" -a"
    print(cov_line, flush=True)
    cov = subprocess.run(cov_line, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

    scov = [line.split("\t") for line in cov.split("\n") if line != '']
    doc = np.mean([int(l[2]) for l in scov])
    boc = len([s for s in scov if int(s[2]) >= 5])/(nv.end-nv.start+1)
    print(f"{scheme}\t{id}\t{doc}\t{boc}", file=fout, flush=True)

    if doc >= 10 and boc >= 0.90:

        out_fq = f"{prefix}.{id}.fq"
        bam2fq_line = f"samtools view {bam_path} -h \"{nv.chr}:{nv.start}-{nv.end}\" | samtools consensus - -d 10"
        consensus_seq = subprocess.run(bam2fq_line, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

        consensus_seq = "".join(consensus_seq.split("\n")[1:])
        return consensus_seq

    else:
        cov_warnings.append(id)
        print(f"{prefix}.{id} not sufficiently covered: {doc}\t{boc}\n", flush=True)
        return None


def run_pipeline(args):
    ref_fa = "/home/amorris/BioInf/MLVA_getter/cryptosporidium_parvum.fasta"
    ref_gff = "/home/amorris/BioInf/MLVA_getter/cryptosporidium_parvum.gff"

    var_profiles = read_yaml(args.yaml)

    var_concats = defaultdict(list)

    bams_list = args.bams

    scheme_loci_ids = {}
    for scheme, s_data in var_profiles.items():
        scheme_loci = read_bed(s_data["bed"])
        scheme_loci_ids[scheme] = scheme_loci

    #     for id, nv in scheme_loci.items():
    #         get_locus_line = f"echo \"{nv.chr}\t{nv.start}\t{nv.end}\" | bedtools getfasta -fi {ref_fa} -bed - > {id}.fa"
    #         subprocess.call(get_locus_line, shell=True)

    loci_dict = defaultdict(lambda: defaultdict(dict))

    for bam_path in bams_list:
        bam_prefix = os.path.basename(bam_path).split("_grouped.bam")[0]
        print(bam_path, bam_prefix)

        with open(f"{bam_prefix}_ST_stats.csv", 'w') as fout:
            print("scheme\tlocus\tDOC\tBOC", file=fout, flush=True)

            for scheme, s_data in var_profiles.items():

                scheme_loci = read_bed(s_data["bed"])

                for id, nv in scheme_loci.items():
                    consensus_fasta = get_locus_seq(f"./{id}.fa", scheme, id, nv, bam_path, bam_prefix, fout)

                    if consensus_fasta == None:
                        print(f"gene extraction for {scheme} : {id} : {bam_prefix} failed.\n", flush=True)
                        continue

                    loci_dict[scheme][bam_prefix][id] = consensus_fasta

    for scheme, sample_cc in loci_dict.items():
        numlocs = len(scheme_loci_ids[scheme])

        concat_fasta = f"{scheme}_concats.fasta"
        var_concats[scheme] = concat_fasta

        with open(concat_fasta, "w") as fout:
            for sample_id, loci in sample_cc.items():
                n = len(loci)

                if n == numlocs:
                    ## conserve sequence order
                    concat_seq = ""
                    for nv in scheme_loci_ids[scheme]:
                        concat_seq += loci[nv]

                    print(f">{sample_id}\n{concat_seq}\n", file=fout, flush=True)
                else:
                    print(f"{sample_id} only covered {n}/{numlocs} of the target loci in scheme {scheme} and will not be included in the alignment.", flush=True)

    # ## run align and build tree
    aln_and_tree(var_concats)

    make_JSON(var_profiles)


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('yaml', action='store', help='YAML file with variant profile definitions.')

    # parser.add_argument('--fastq', action='store', nargs=2, help='FASTQ files to map.')
    parser.add_argument('--bams', action='store', nargs="*", default=None, help='input BAM files.')

    args = parser.parse_args(argv)

    return args

if __name__=="__main__":
    args = parse_args(sys.argv)
    run_pipeline(args)
