// modules for mapping and processing fastq files

process checkFqValidity {
  /**
  * Check if fastq files are valid
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing", mode: 'copy', overwrite: 'true', pattern: '*.err'

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
  printf OK
  touch ${error_log}
  """
}

process fastQC {
  /**
  * run FastQC on each fastq pair and copy reports to a subdir in outputdir
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/fastQC", mode: 'copy'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)

  output:
  path("*", emit: fastQC_report)

  script:
  fqc_log = "${sample_name}fqc.log"

  """
  cat $fq1 $fq2 > ${sample_name}.fq.gz
  fastqc ${sample_name}.fq.gz 1> $fqc_log
  rm ${sample_name}.fq.gz
  """

  stub:
  report="${sample_name}.html"

  """
  touch ${report}
  """
}

process multiQC {
  /**
  * TODO: this doesn't work at all, I don't know why
  */

  publishDir "${params.output_dir}/QC_reports", mode: 'copy'

  memory '5 GB'

  input:
  path(fastQC_reports)

  output:
  path("fastqc_paths.txt", emit: multiQC_report)

  script:
  """
  echo ${fastQC_reports} >> fastqc_paths.txt
  """
}

process trimGalore {
  /**
  * run trimGalore on fastq pair
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing/Trim", mode: 'copy'

  memory '10 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)

  output:
  tuple path("${fq1}_trimming_report.txt"), path("${fq2}_trimming_report.txt"), emit: tg_reports
  tuple path("*_val_1.*"), path("*_val_2.*"), emit: tg_fqs

  script:
  tg_log = "${sample_name}_tg.log"
  """
  trim_galore --paired $fq1 $fq2 1> $tg_log
  """

  stub:
  rep1 = "${fq1}_trimming_report.txt"
  rep2 = "${fq2}_trimming_report.txt"
  trim_fq1 = "${sample_name}_val_1.fq.gz"
  trim_fq2 = "${sample_name}_val_2.fq.gz"

  """
  touch ${rep1}
  touch ${rep2}
  touch ${trim_fq1}
  touch ${trim_fq2}
  """
}

process mapBT2 {
  /**
  * map reads to a reference genome using Bowtie2
  */

  tag { sample_name }

  cpus 8
  memory '15 GB'

  publishDir "${params.output_dir}/$sample_name/preprocessing/Map", mode: 'copy'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(trimmed_fq1), path(trimmed_fq2)
  path(ref_bt2index)

  output:
  path("${sample_name}.sorted.bam"), emit: bam
  path("${sample_name}_alnStats.txt"), emit: map_stats


  script:
  bam = "${sample_name}.sorted.bam"
  bai = "${sample_name}.bam.bai"
  stats = "${sample_name}_alnStats.txt"

  """
  echo $ref_bt2index
  bowtie2 --very-sensitive -p ${task.cpus} -x ${workflow.launchDir}/${params.output_dir}/REFDATA/${params.genome}/${params.genome} -1 $fq1 -2 $fq2 2> ${sample_name}_alnStats.txt | samtools view -f 4 -Shb - | samtools sort - -o ${bam}
  """

  stub:
  bam = "${sample_name}.sorted.bam"
  stats_txt = "${sample_name}_alnStats.txt"

  """
  touch ${bam}
  touch ${stats_txt}
  """
}

process deduplication {
  /**
  * run picard to perform deduplication on reads
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing/Dedup", mode: 'copy'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(bam)

  output:
  path("${sample_name}_dedup.bam"), emit: dedup_bam
  path("${sample_name}_metrics.txt"), emit: dedup_metrics

  script:
  dedup_log = "${sample_name}_dedup.log"
  """
  java -jar /usr/local/bin/picard.jar MarkDuplicates INPUT=${bam} OUTPUT=${sample_name}_dedup.bam METRICS_FILE=${sample_name}_metrics.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 1> $dedup_log
  """

  stub:
  dedup_bam = "${sample_name}_dedup.bam"
  metrics = "${sample_name}_metrics.txt"

  """
  touch ${dedup_bam}
  touch ${metrics}
  """
}

process qualimap {
  /**
  * TODO: this doesn't seem to behave as expected. Claims no reads map.
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing/qualimap", mode: 'copy'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(dedup_bam)

  output:

  script:
  """
  qualimap bamqc -bam ${dedup_bam} -outdir ./ -outformat HTML
  """
}
