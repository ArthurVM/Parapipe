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


def make_multifasta(scheme_id, concats):

    out_path = f"{scheme_id}.concats.fasta"
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


def aln_and_tree(scheme_id, concat_fasta):
    """ Requires clustalw and iqtree2 in PATH
    """

    # for scheme_id, multifasta in var_concats.items():
    # out_path = make_multifasta(scheme_id, concats)

    aln_line = f"clustalw2 -INFILE={concat_fasta} -ALIGN -OUTPUT=nexus -TYPE=DNA -OUTFILE={scheme_id}.nxs"
    subprocess.call(aln_line, shell=True)

    tree_line = f"iqtree2 -s {scheme_id}.nxs -b 100"
    subprocess.call(tree_line, shell=True)


def make_JSON(var_profiles, scheme_dict):

    for scheme_id, bed in var_profiles.items():
        scheme_dict[scheme_id]["bed"] = bed

        contree = scheme_id + ".nxs.contree"
        iqtree_log = scheme_id + ".nxs.iqtree"
        mldist = scheme_id + ".nxs.mldist"
        nxs = scheme_id + ".nxs"

        with open(contree, 'r') as fin:
            data = fin.read().strip("\n")
            scheme_dict[scheme_id]["contree"] = data

        with open(mldist, 'r') as fin:
            mlbox = []
            lines = fin.readlines()

            ## deal with failed scheme analyses
            if lines[0].strip("\n") == "None":
                scheme_dict[scheme_id]["mldist"] = "None"

            else:
                for line in lines:
                    sline = [l for l in line.strip("\n").split(" ") if l!='']
                    mlbox.append(sline)

                scheme_dict[scheme_id]["mldist"] = mlbox

    with open("typing_report.json", 'w') as fout:
        json.dump(scheme_dict, fout, indent=4)


def get_locus_seq(args, ref, scheme, id, nv, bam_path, prefix):

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

    if doc >= args.doc and boc >= args.boc:

        out_fq = f"{prefix}.{id}.fq"
        bam2fq_line = f"samtools view {bam_path} -h \"{nv.chr}:{nv.start}-{nv.end}\" | samtools consensus - -d 10"
        consensus_seq = subprocess.run(bam2fq_line, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

        consensus_seq = "".join(consensus_seq.split("\n")[1:])
        return consensus_seq, doc, boc, cov_warnings

    else:
        cov_warnings.append(id)
        print(f"{prefix}.{id} not sufficiently covered: {doc}\t{boc}\n", flush=True)
        return None, doc, boc, cov_warnings


def make_empty_out(scheme_id):
    outfiles = [f"{scheme_id}.nxs.contree", f"{scheme_id}.nxs.iqtree", f"{scheme_id}.nxs.mldist", f"{scheme_id}.nxs", "snp_tree.nwk"]

    for outfile in outfiles:
        with open(outfile, 'w') as fout:
            print(f"None", file=fout)


def run_pipeline(args):
    ref_fa = "/home/amorris/BioInf/MLVA_getter/cryptosporidium_parvum.fasta"
    ref_gff = "/home/amorris/BioInf/MLVA_getter/cryptosporidium_parvum.gff"

    var_profiles = read_yaml(args.yaml)
    var_concats = defaultdict(list)
    bams_list = args.bams

    scheme_dict = defaultdict(lambda: defaultdict(dict))
    scheme_loci_ids = {scheme: read_bed(s_data["bed"]) for scheme, s_data in var_profiles.items()}
    loci_dict = defaultdict(lambda: defaultdict(dict))

    for bam_path in bams_list:
        bam_prefix = os.path.basename(bam_path).split("_grouped.bam")[0]

        for scheme, scheme_loci in scheme_loci_ids.items():
            for id, nv in scheme_loci.items():

                print(f"Processing scheme: {scheme}, id: {id}, bam_prefix: {bam_prefix}")

                ## initialise dictionary structure  
                if bam_prefix not in scheme_dict[scheme]["locus_coverage"]:
                    scheme_dict[scheme]["locus_coverage"][bam_prefix] = {}
                if id not in scheme_dict[scheme]["locus_coverage"][bam_prefix]:
                    scheme_dict[scheme]["locus_coverage"][bam_prefix][id] = {}

                consensus_fasta, doc, boc, cov_warnings = get_locus_seq(args, f"./{id}.fa", scheme, id, nv, bam_path, bam_prefix)

                scheme_dict[scheme]["locus_coverage"][bam_prefix][id]["BOC"] = boc
                scheme_dict[scheme]["locus_coverage"][bam_prefix][id]["DOC"] = doc
                scheme_dict[scheme]["locus_coverage"][bam_prefix][id]["warnings"] = cov_warnings

                if consensus_fasta is None:
                    print(f"gene extraction for {scheme} : {id} : {bam_prefix} failed.\n", flush=True)
                    continue

                loci_dict[scheme][bam_prefix][id] = consensus_fasta

    ## check that there is enough sequence data to proceed
    ## there may not be if all samples were low quality and therefore no loci for any schemes were covered
    if len(loci_dict) == 0:
        print(f"All schemes failed due to lack of coverage of constituent loci across the dataset.")
        for scheme in scheme_loci_ids:
            make_empty_out(scheme)
    
    else:
        ## check if any schemes failed
        for scheme in scheme_loci_ids:
            if scheme not in loci_dict:
                print(f"Scheme {scheme} failed due to lack of coverage of constituent loci across the dataset.")
                make_empty_out(scheme)

        ## proceed with successful schemes
        for scheme, sample_cc in loci_dict.items():
            numlocs = len(scheme_loci_ids[scheme])
            concat_fasta = f"{scheme}_concats.fasta"
            var_concats[scheme] = concat_fasta
            seqs_in_concat_fasta = 0

            with open(concat_fasta, "w") as fout:
                for sample_id, loci in sample_cc.items():
                    n = len(loci)

                    if n == numlocs:
                        concat_seq = "".join(loci[nv] for nv in scheme_loci_ids[scheme] if loci[nv] and len(loci[nv]) >= 10)
                        print(f">{sample_id}\n{concat_seq}\n", file=fout, flush=True)
                        seqs_in_concat_fasta += 1
                    else:
                        print(f"{sample_id} only covered {n}/{numlocs} of the target loci in scheme {scheme} and will not be included in the alignment.", flush=True)

            if seqs_in_concat_fasta >= 3:
                aln_and_tree(scheme, concat_fasta)
            else:
                print(f"Fewer than 3 ({seqs_in_concat_fasta}) sequences detected for scheme {scheme}. Cannot run alignment or tree building.")
                make_empty_out(scheme)

    make_JSON(var_profiles, scheme_dict)


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('yaml', action='store', help='YAML file with variant profile definitions.')

    # parser.add_argument('--fastq', action='store', nargs=2, help='FASTQ files to map.')
    parser.add_argument('--bams', action='store', nargs="*", default=None, help='input BAM files.')
    parser.add_argument('--doc', action='store', default=5, type=int, help='Minimum depth of coverage across a locus to include in scheme. Default=5.')
    parser.add_argument('--boc', action='store', default=0.9, type=float, help='Minimum breadth of coverage across a locus to include in scheme. Default=0.9.')

    args = parser.parse_args(argv)

    return args

if __name__=="__main__":
    args = parse_args(sys.argv)
    run_pipeline(args)
