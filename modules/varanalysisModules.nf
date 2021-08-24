// modules for identifying variation in mapped reads

process runSNP {
  /**
  * Call snps in mapped reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(dedup_bam)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("${sample_name}.var.vcf"), emit: vcf

  script:
  error_log = "${sample_name}.err"

  """
  samtools mpileup --skip-indels --BCF -f ${fasta} ${dedup_bam} | bcftools call --skip-variants indels -m -O v --variants-only -o ${sample_name}.var.vcf -
  """

  stub:
  error_log  = "${sample_name}.err"

  """
  touch ${sample_name}.var.vcf
  """
}

process preprocessForPlotting {
  /**
  * preprocesses FASTA and VCF files for plotting
  */

  tag { sample_name }

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("tmp.descStripped.fasta"), emit: preprocessed_fasta
  path("vcf_split"), emit: vcf_path

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  echo "Preprocessing ${fasta}..."
  python3 ${scripts}/removeFASTAdesc.py ${fasta}
  echo "Splitting ${vcf} by chromosome ID..."
  python3 ${scripts}/splitVCF.py ${vcf}
  """

  stub:
  """
  touch tmp_descStripped.fasta
  mkdir vcf_split
  """
}

process plotSNP {
  /**
  * Plot SNPs across each chromosome
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf_path)
  path(preprocessed_fasta)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("SNPDist.pdf"), emit: pdf

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  echo "Creating plots..."
  Rscript ${scripts}/SNPdist.R -v ${vcf_path} -f ${preprocessed_fasta} -g ${gff}
  """

  stub:
  """
  touch SNPDist.pdf
  """
}

process findSTRs {
  /**
  * Find Short Tandem Repeats using TRF
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/STR", mode: 'copy', overwrite: 'true', pattern: "*.dat"

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(pilon_fasta)

  output:
  path("trf.log")

  script:
  trf_log = "trf.log"
  """
  trf4.10.0-rc.2.linux64.exe ${pilon_fasta} 2 7 7 80 10 50 500 -d -h &> ${trf_log}
  """

  stub:

  """
  """
}
