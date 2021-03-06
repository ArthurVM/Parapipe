// modules for mapping and processing fastq files

process checkFqValidity {
  /**
  * Check if fastq files are valid
  */

  tag { sample_name }

  errorStrategy 'ignore'

  publishDir "${params.output_dir}/$sample_name/preprocessing", mode: 'copy', overwrite: 'true', pattern: '*.err'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)

  output:
  tuple val(sample_name), path(fq1), path(fq2), emit: checkValidity_fqs
  path("${sample_name}.err", emit: checkValidity_log)

  script:
  error_log = "${sample_name}.err"

  """
  is_ok=\$(fqtools validate $fq1 $fq2)

  if [ \$is_ok == 'OK' ]; then printf 'OK' && printf "" >> ${error_log}; else echo "error: sample did not pass fqtools validation check\n \$is_ok" >> ${error_log}; fi
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

  errorStrategy 'ignore'

  publishDir "${params.output_dir}/${sample_name}/preprocessing", mode: 'copy'

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

  publishDir "${params.output_dir}", mode: 'copy'

  memory '5 GB'

  input:
  path(fastqc_files)

  output:
  path("multiqc_report.html", emit: multiQC_report)

  script:
  """
  echo ${fastqc_files}
  multiqc ./
  """

  stub:
  multiQC_report="multiqc_report.html"
  """
  touch ${multiQC_report}
  """
}

process trimGalore {
  /**
  * run trimGalore on fastq pair
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing/Trim", mode: 'copy', pattern: '*.trimming_report.txt'

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

process map2Ref {
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
  bowtie2 --very-sensitive -p ${task.cpus} -x ${workflow.launchDir}/${params.output_dir}/REFDATA/${params.genome}/${params.genome} -1 $fq1 -2 $fq2 2> ${sample_name}_alnStats.txt | samtools view -h - | samtools sort - -o ${bam}
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

process gini {
  /**
  * calculate the gini of mapped read coverage
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/preprocessing/gini", mode: 'copy'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(bam)

  output:
  path "${sample_name}.GG", emit: GG
  path "${sample_name}.doc.bed", emit: doc_bed

  script:
  """
  samtools depth -a ${bam} > ${sample_name}.doc.bed
  python3 /NGS-Gini-Analysis-Toolkit-0.1-alpha/src/gini.py ${sample_name}.doc.bed -G 5 1000 50 > ${sample_name}.GG
  """

  stub:
  gg_file = "${sample_name}.GG"
  """
  touch ${gg_file}
  """
}

process summarise {
  /**
  * calculate the gini of mapped read coverage
  */
  tag { sample_name }

  publishDir "${params.output_dir}/mapping_stats", mode: 'copy'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(bam)
  path(doc_bed)
  path(GG)

  output:
  path "${sample_name}_mapstats.json", emit: mapstats_json

  script:
  """
  samtools index ${bam}
  samtools stats ${bam} > tmp.stats
  python3 ${baseDir}/scripts/parse_samtools_stats.py ${bam} tmp.stats ${doc_bed} ${GG} > ${sample_name}_mapstats.json
  """

  stub:
  mapstats_json = "${sample_name}_mapstats.json"
  """
  touch ${mapstats_json}
  """
}
