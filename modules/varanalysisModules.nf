// modules for identifying variation in mapped reads

process runSNP {
  /**
  * Call snps in mapped reads
  * TODO: something weird is happening and it is not reporting any SNPs, I dont know why this is
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantCall", mode: 'copy', overwrite: 'true'

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
