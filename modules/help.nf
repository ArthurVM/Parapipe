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


  Profiles:
  singularity           run with singularity image
  """.stripIndent()
}