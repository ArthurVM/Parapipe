// modules for identifying target loci in unmapped reads

process runBlooMine {
  /**
  * Process STRs discovered by Tandem Repeats Finder using SNPEff
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/MOI", mode: 'copy', overwrite: 'true'

  cpus 2
  memory '25 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(trimmed_fq1), path(trimmed_fq2)
  path(target_flanks)

  output:
  path("${sample_name}*.txt"), emit: moi_report

  script:
  """
  zcat ${trimmed_fq1} ${trimmed_fq2} > ${sample_name}_UNP.fq
  python3 /BlooMine/BlooMine.py ${target_flanks} ${sample_name}_UNP.fq -p ${sample_name} -t 1
  """

  stub:
  ann_vcf = "${sample_name}.target1_BMfiltered.report.txt"
  """
  touch ${ann_vcf}
  """
}
