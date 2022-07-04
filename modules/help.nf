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
    --genome            the ID of the reference genome to download and map against.
    --ref_scaffold      perform reference driven scaffolding using ABACAS. NOTE: It is unadvised to use this if the reference is not arranged in complete chromosomes.
    --only_SNP          only run SNP calling modules and exit before assembly and downstream analysis.
    --target_flanks     flaking regions for target loci to extract for the purpose of phylogenetic analysis.

  Profiles:
  singularity           run with singularity image
  """.stripIndent()
}
