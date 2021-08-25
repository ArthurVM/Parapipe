// modules for identifying STRs in mapped reads

process findSTRs {
  /**
  * Find Short Tandem Repeats using TRF
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/STR", mode: 'copy', overwrite: 'true'

  memory '12 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(pilon_fasta)

  output:
  tuple path("trf.stdout"), path("trf.stderr"), emit: trf_log
  path("${pilon_fasta}.2.7.7.80.10.50.500.dat"), emit: trf_dat

  script:
  trf_stdout = "trf.stdout"
  trf_stderr = "trf.stderr"
  // There is some weirdness going on with TRF, where it is successfully running but producing non 0 exit status
  // trap is used here as a temporary fix
  """
  trap 'if [[ \$? != 0 ]]; then echo OVERRIDE EXIT; exit 0; fi' EXIT
  trf4.10.0-rc.2.linux64.exe ${pilon_fasta} 2 7 7 80 10 50 500 -d -h 1> ${trf_stdout} 2> ${trf_stderr}
  """

  stub:
  trf_dat = "${pilon_fasta}.2.7.7.80.10.50.500.dat"
  """
  touch ${trf_dat}
  """
}
