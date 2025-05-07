<p align="center">
    <img src="resources/PP_logo.png" alt="Parapipe logo" width="500">
</p>

# Parapipe: A Bioinformatics Pipeline for Genomic Analysis of Eukaryotic Pathogens from NGS Data

## Overview

Parapipe is a modular, ISO-accreditable bioinformatics pipeline designed to process Eukaryotic Pathogen next-generation sequencing (NGS) datasets. It provides a robust, scalable, and reproducible framework for genomic surveillance, outbreak tracking, and epidemiological investigations.

Parapipe supports paired end Next-Generation Sequencing data and performs thorough Quality-Control, ensuring high-quality genomic insights even when full genome coverage is not achieved.

## Features

- **Quality Control**: Automated read pre-processing, filtering, and QC reporting ([fastp](https://github.com/OpenGene/fastp), [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/), [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)).
- **Variant Calling & SNP Analysis**: Variant calling and high-resolution SNP space analysis ([Samtools](http://www.htslib.org/), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [FreeBayes](https://github.com/ekg/freebayes)).
- **Multiplicity of Infection (MOI) Characterisation**: Identifies mixed infections and within-host diversity ([PyMOI](https://github.com/ArthurVM/PyMOI), [Moimix](https://github.com/bahlolab/moimix)).
- **Phylogenomic Analysis**: Generates SNP-based phylogenetic trees and clustering (Ardal, [IQtree2](https://github.com/iqtree/iqtree2)).
- **Modular and Scalable**: Built using [**Nextflow DSL2**](https://www.nextflow.io/), containerized with [**Singularity**](https://sylabs.io/guides/3.0/user-guide/), ensuring portability and reproducibility.

## Installation

### Prerequisites
Bioinformatic tools executed within this pipeline itself are containerised using singularity/apptainer, but Go, Nextflow, and Singularity/Apptainer themselves will need to be installed.
- [Nextflow](https://www.nextflow.io/)
- [Singularity](https://sylabs.io/guides/3.0/user-guide/)
- [Go](https://go.dev/)

### Installation Steps

1. Clone the repository:
   ```
   git clone https://github.com/ArthurVM/Parapipe.git
   ```
2. Build Singularity containers:
   ```
   cd Parapipe/singularity
   bash singularity_build.sh
   cd ..
   ```

## Run Parapipe
Parapipe runs using the standard Nextflow command line syntax. There are a number of parameters which can be assigned and adjusted using command line arguments, but you can run Parapipe out of the box.

### Quick Start
For those who want to just run Parapipe out of the box with default parameters, just feed it a directory containing paired end fastqs, a suffix for how fastqs can be paired, a reference species to map to, and an output directory prefix. 

Lets say we have a dataset of paired end Cryptosporidium parvum NGS files in our directory, `/user/my_cparvum_data/`, which looks like:
```
/user/my_cparvum_data/
            ├── sample1_1.fastq.gz
            ├── sample1_2.fastq.gz
            ├── sample2_1.fastq.gz
            ├── sample2_2.fastq.gz
            ├── sample3_1.fastq.gz
            ├── sample3_2.fastq.gz
            ├── sample4_1.fastq.gz
            └── sample4_2.fastq.gz
```
The `--pattern` argument refers to the suffix with which we can pair files. In this case the suffix for forward reads are `_1.fastq.gz`, and reverse reads `_2.fastq.gz`. Consequently, we can use bash syntax to pair these files by setting `--pattern` to `'*_{1,2}.fastq.gz'`, which captures any files which end in `_1.fastq.gz` or `_2.fastq.gz` using the `*` wildcard, and pairs them by their shared prefix (e.g. `sample1`). This prefix will then be used as the ID for this sample during the run.

In the Parapipe directory, simply run:
```
nextflow run ./parapipe.nf --input_dir /user/my_cparvum_data/ --pattern '*_{1,2}.fastq.gz' --ref cryptosporidium_parvum --output_dir ./my_cparvum_run -profile singularity
```

### Advanced Usage
There are a small number of parameters which can be tuned, determining which samples are included in phylogenomic analysis, and how alleles are filtered.


```
❯ nextflow run ./parapipe.nf --help

Usage:
  nextflow run main.nf --input_dir [fullPath] --pattern [regexPattern] --genome [refGenomeID]

Description:
  General purpose parasite genomics pipeline. DEVELOPMENT VERSION.

Mandatory Arguments:
  --input_dir         the full path to the directory containing raw read files in paired-end
                      fastq format.
  --output_dir        output directory. Module output files will be written to subdirectories
                      within this parent directory.
  --pattern           regex pattern to match pairs of fastq files.
  --ref               the species ID of the reference genome to download and map against.
  --yaml              YAML file containing typing profiles.

Optional Arguments:
  --database          the path to a directory containing vcf files for constructing a
                      phylogenetic tree.
  --mincov            the minimum fraction of the genome which must be covered to a depth of
                      5x to include a sample in phylogenetic analysis. Default=0.8.
  --missing           the maximum allele missingness to include SNPs in MOI analysis. Default=0.1.
  --maf               the minimum minor allele frequency to include SNPs in MOI analysis. Default=0.05.
  --mac               the minimum minor allele depth to include SNPs in MOI analysis. Default=5.

Profiles:
singularity           run with singularity image
```

## Output
Parapipe will output a directory which contains key intermediate and analysis files which can be used for further analysis. A report is produced and deposited in the output directory as `Parapipe_report.html`.

## Suggestions and Requests
If you have any suggestions for functionality, requests for directions on usage or report interpretation, or if you would just like to reach out to discuss Parapipe then please don't hesitate to contact me at morrisa28@cardiff.ac.uk. Alternatively you may wish to post an issue using the GitHub issues tab.

## Citing Parapipe
If you use Parapipe to perform analysis of your dataset, please cite the preprint at https://www.biorxiv.org/content/10.1101/2025.01.16.633364v1
