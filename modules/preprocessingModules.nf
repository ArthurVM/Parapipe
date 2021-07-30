// modules for mapping and processing fastq files

process checkFqValidity {
    /**
    * @QCcheckpoint confirm that fqtools validates both fastqs
    */

    tag { sample_name }

    publishDir "${params.output_dir}/$sample_name", mode: 'copy', overwrite: 'true', pattern: '*.err'

    memory '5 GB'

    input:
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path("${sample_name}.err", emit: checkValidity_log)

    script:
    error_log = "${sample_name}.err"

    """
    is_ok=\$(fqtools validate $fq1 $fq2)

    if [ \$is_ok == 'OK' ]; then printf 'OK' && printf "" >> ${error_log}; else echo "error: sample did not pass fqtools validation check" >> ${error_log}; fi
    """

    stub:
    error_log  = "${sample_name}.err"

    """
    printf ${params.checkFqValidity_isok}
    touch ${error_log}
    """
}

process map2Ref {
  /**
  * map reads in fastq format to a reference genome formatted as a BT2 index set
  */

  cpus 8
  memory '15 GB'

  publishDir:

  input:
  val(saple_name)
  path(fqDir)
  path(ref_bt2index)

  output:
  path("${sample_name}.bam")
  path("${sample_name}.bai")
  path("${sample_name}_alnStats.txt")

  script:
  bam = "${sample_name}.bam"
  bai = "${sample_name}.bam.bai"
  stats = "${sample_name}_alnStats.txt"

  """

  """
}
