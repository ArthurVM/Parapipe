import os
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


def get_locseqs(fasta, locdict):

    seqdict = defaultdict(str)

    recs = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    for id, loc in locdict.items():
        seq = str(recs[loc.chr].seq[loc.start-1:loc.end])
        seqdict[id] = seq

    with open("targetLoci.fasta", 'w') as fout:
        for id, rec in seqdict.items():
            fout.write(f">{id}\n{seq}\n")

    return seqdict, "./targetLoci.fasta"


def run_bt2(ref, fq1, fq2):

    prefix = os.path.basename(fq1).split("_1")[0]
    ref_prefix = os.path.basename(ref).split(".fasta")[0]

    if not os.path.exists(ref_prefix):
        os.mkdir(ref_prefix)

    build_runline = f"bowtie2-build {ref} {ref_prefix}/{ref_prefix} > buildlog.txt"
    print(build_runline)
    subprocess.call(build_runline, shell=True)

    map_runline = f"bowtie2 --very-sensitive -p 6 -x {ref_prefix}/{ref_prefix} -1 {fq1} -2 {fq2} 2> {prefix}_alnStats.txt | samtools view -h - | samtools sort - -o {prefix}.sorted.bam"
    print(map_runline)
    subprocess.call(map_runline, shell=True)

    return f"{prefix}.sorted.bam"


def make_var_fasta(vcf, fasta, locdict, prefix, profile_id, qual_t):

    var_sequences = defaultdict(str)
    ref_sequences = defaultdict(str)

    if not os.path.exists(f"{vcf}.tbi"):
        vcf_index_line = f"bcftools index {vcf}"
        print(vcf_index_line)
        subprocess.call(vcf_index_line, shell=True)

    vcf = pysam.VariantFile(vcf)
    ref_recs = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    ## Process each gene
    for id, nv in locdict.items():
        ## Extract the gene sequence from the reference genome
        gene_sequence = ref_recs[nv.chr][nv.start - 1:nv.end]
        var_sequence = ref_recs[nv.chr][nv.start - 1:nv.end]

        ## initialise a buffer to adjust variant positions for previous sequence introductions
        buffer = 0

        ## Iterate through variants in the VCF file that overlap with the gene
        for record in vcf.fetch(nv.chr, nv.start, nv.end):

            ## ensure the variant exceeds the given quality threshold
            if record.qual >= qual_t:

                ## Check if the variant is within the gene
                if nv.start <= record.pos <= nv.end:
                    ## Get the alternate allele (variant) sequence
                    ref_allele = record.alleles[0]
                    alt_allele = record.alleles[1]

                    ## Calculate the position of the variant within the gene
                    variant_pos = record.pos - nv.start + buffer

                    ## Replace the reference allele with the alternate allele
                    ref_locus_sequence = var_sequence.seq[variant_pos: variant_pos + len(ref_allele)]

                    # print(ref_locus_sequence, record.pos, ref_allele, alt_allele, buffer)

                    if ref_locus_sequence != ref_allele:
                        print("Something went wrong")
                        print(id)
                        print(f"found_R={ref_locus_sequence} R={ref_allele} A={alt_allele} buffer={buffer} pos={record.chrom}:{record.pos}")
                        sys.exit(1)

                    # print(ref_locus_sequence, ref_allele, alt_allele, buffer)
                    var_sequence = var_sequence[:variant_pos] + alt_allele + var_sequence[variant_pos + len(ref_allele):]
                    buffer += len(alt_allele)-len(ref_allele)

        ## Store the modified gene sequence
        var_sequences[id] = var_sequence
        ref_sequences[id] = gene_sequence

    ## write sequence files
    with open(f"{prefix}.{profile_id}.fasta", 'w') as fout, \
    open(f"{prefix}.{profile_id}_concat.fasta", 'w') as cc_fout, \
    open(f"REF_{profile_id}.concat.fasta", 'w') as cc_ref_fout:
        print(f">{prefix}", file=cc_fout)
        print(f">REF_{profile_id}", file=cc_ref_fout)
        for id, rec in var_sequences.items():
            print(f">{id}\n{rec.seq}", file=fout)
            print(rec.seq, end="", file=cc_fout)
            print(ref_sequences[id].seq, end="", file=cc_ref_fout)

    return f"{prefix}.{profile_id}_concat.fasta"


def freebayes(fasta, bam, bam_prefix):
    out_vcf = f"./{bam_prefix}.vcf"

    if not os.path.exists(out_vcf + ".gz"):
        ## sed 1: modify sample ID from unknown
        ## sed 2: modify GQ type to string to avoid bcftools type error when encountering a . in this field
        vcall_runline = f"freebayes -p 2 -P 0 -C 2 -F 0.05 --min-coverage 5 --min-repeat-entropy 1.0 -q 13 -m 1 --gvcf --gvcf-dont-use-chunk true --strict-vcf -f {fasta} {bam} > {bam_prefix}.tmp.vcf \
            && sed -i 's/##FORMAT=<ID=GQ,Number=1,Type=Integer/##FORMAT=<ID=GQ,Number=1,Type=String/' {bam_prefix}.tmp.vcf \
            && sed -i 's/unknown/{bam_prefix}/' {bam_prefix}.tmp.vcf"

        print(vcall_runline)
        subprocess.call(vcall_runline, shell=True)

        sort_runline = f"bcftools sort {bam_prefix}.tmp.vcf -o {out_vcf} \
            && bgzip -c ./{out_vcf} > ./{out_vcf}.gz \
            && tabix -f -p vcf ./{out_vcf}.gz"

        print(sort_runline)
        subprocess.call(sort_runline, shell=True)

        ## clean up unnecessary files
        os.remove(out_vcf)
        os.remove(f"{bam_prefix}.tmp.vcf")

    else:
        print(f"Found {out_vcf}.gz")

    return out_vcf + ".gz"


def aln_and_tree(var_concats, t):
    """ Requires clustalw and iqtree2 in PATH
    """

    ## store schemes for which trees were successfully built
    tree_pass = {}

    for profile_id, concats in var_concats.items():
        if len(concats) < 4:
            print(f"Fewer than 4 sequences detected for scheme {profile_id}. Cannot construct tree.")
            tree_pass[profile_id] = False
            continue

        out_path = make_multifasta(profile_id, concats)

        aln_line = f"clustalw2 -INFILE={out_path} -ALIGN -OUTPUT=nexus -TYPE=DNA -OUTFILE={profile_id}.nxs"
        subprocess.call(aln_line, shell=True)

        tree_line = f"iqtree2 -s {profile_id}.nxs -nt {t} -b 1000"
        subprocess.call(tree_line, shell=True)

        tree_pass[profile_id] = True

    return tree_pass


def make_JSON(tree_pass, var_profiles):

    json_dict = defaultdict(dict)

    for profile_id, bed in var_profiles.items():

        tree_bool = tree_pass[profile_id]
        
        ## populate the json if there were too few sequences to build a tree
        if not tree_bool:
            json_dict[profile_id]["bed"] = "NA"
            json_dict[profile_id]["contree"] = "NA"
            json_dict[profile_id]["mldist"] = "NA"
            continue

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


def get_loc_cov(scheme, locdict, bam_path, prefix, fout):
    """ get the coverage over each region within the typing scheme
    """

    cov_warnings = []

    for locus_id, nv in locdict.items():

        if not os.path.exists(f"{bam_path}.bai"):
            bam_index_line = f"samtools index {bam_path}"
            print(bam_index_line)
            subprocess.call(bam_index_line, shell=True)

        cov_line = f"samtools depth {bam_path} -r \"{nv.chr}:{nv.start}-{nv.end}\" -a"
        print(cov_line)
        cov = subprocess.run(cov_line, shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')

        scov = [line.split("\t") for line in cov.split("\n") if line != '']
        doc = np.mean([int(l[2]) for l in scov])
        boc = len([s for s in scov if int(s[2]) >= 5])/(nv.end-nv.start+1)
        print(f"{scheme}\t{locus_id}\t{doc}\t{boc}", file=fout)

        if doc <= 5 or boc <= 0.75:
            print(f"WARNING: coverage over specified region is low.")
            print(f"Scheme:\t{scheme}\nLocus:\t{locus_id}\tDOC:\t{doc}\tBOC:\t{boc}\n")
            cov_warnings.append(locus_id)

    lowq_loci = '\t'.join(cov_warnings)
    print(f"WARNINGS\t{scheme}\t{lowq_loci}", file=fout)

    return cov_warnings


def run_pipeline(args):

    bam_suffix = "_grouped.bam"
    var_profiles = read_yaml(args.yaml)

    var_concats = defaultdict(list)

    bams_list = [bam for bam in os.listdir(args.bams) if bam.endswith(bam_suffix)]

    if len(bams_list) == 0:
        print(f"No bams found in {args.bams} of form *{bam_suffix}")
        sys.exit(666)

    for bam in bams_list:
        bam_prefix = os.path.basename(bam).split(bam_suffix)[0]
        bam_path = os.path.join(args.bams, bam)

        with open(f"{bam_prefix}_ST_stats.csv", 'w') as fout:
            print("scheme\tlocus\tDOC\tBOC", file=fout)

            vcf = freebayes(args.fasta, bam_path, bam_prefix)

            for scheme, subdict in var_profiles.items():
                bed = subdict["bed"]
                missing = subdict["missing"]
                locdict = read_bed(bed)

                cov_warnings = get_loc_cov(scheme, locdict, bam_path, bam_prefix, fout)

                ## check that the number of low coverage or missing alleles does not exceed the missing threshold provided in the yaml
                if len(cov_warnings) > missing:
                    print(f"ALLELE COV WARNING\n\tscheme : {scheme}\n\tlow cov alleles : {cov_warnings}\n\tnan_threshold : {missing}\n")
                    continue

                concat_fasta = make_var_fasta(vcf, args.fasta, locdict, bam_prefix, scheme, args.qual)
                var_concats[scheme].append(concat_fasta)

    ## run align and build tree
    tree_pass = aln_and_tree(var_concats, args.threads)

    make_JSON(tree_pass, var_profiles)


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('yaml', action='store', help='YAML file with variant profile definitions.')
    parser.add_argument('fasta', action='store', help='FASTA file used as reference to construct the input BAM file.')

    # parser.add_argument('--fastq', action='store', nargs=2, help='FASTQ files to map.')
    parser.add_argument('--bams', action='store', default=None, help='Directory containing input BAM files.')
    parser.add_argument('--db', action='store', default=None, help='Directory containing FreeBayes vcf files to incude in phylogenomic analysis.')
    parser.add_argument('--qual', action='store', default=30, help='Quality score threshold for introducing a variant into a sequence.')
    parser.add_argument('--threads', action='store', default=1, help='The number of threads to use for iqtree2. Default=1')

    args = parser.parse_args(argv)

    return args


if __name__=="__main__":
    args = parse_args(sys.argv)
    run_pipeline(args)
