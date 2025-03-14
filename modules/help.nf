def printHelp() {
log.info"""
  Usage:
    nextflow run main.nf --input_dir [fullPath] --pattern [regexPattern] --genome [refGenomeID]

  Description:
    General purpose parasite genomics pipeline. DEVELOPMENT VERSION.

  Mandatory Arguments:
    --input_dir         the full path to the directory containing raw read files in paired-end fastq format.
    --output_dir        output directory. Module output files will be written to subdirectories within this parent directory.
    --pattern           regex pattern to match pairs of fastq files.
    --ref               the species ID of the reference genome to download and map against.
    --yaml              YAML file containing typing profiles.

  Optional Arguments:
    --read_n_threshold  the minimum number of read pairs in a paired end FASTQ run required for analysis.
    --database          the path to a directory containing vcf files for constructing a phylogenetic tree.
    --mincov            the minimum fraction of the genome which must be covered to a depth of 5x to include a sample in phylogenetic analysis. Default=0.8.
    --missing           the maximum allele missingness to include SNPs in MOI analysis. Default=0.1.
    --maf               the minimum minor allele frequency to include SNPs in MOI analysis. Default=0.05.
    --mac               the minimum minor allele depth to include SNPs in MOI analysis. Default=5.

  Profiles:
  singularity           run with singularity image
  """.stripIndent()
}
