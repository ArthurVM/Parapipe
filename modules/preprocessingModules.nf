// modules for mapping and processing fastq files


process checkFqValidity {
  /**
  * Check if fastq files are valid
  */

  tag { sample_name }

  errorStrategy 'ignore'

  publishDir "${params.output_dir}/${sample_name}/preprocessing", mode: 'copy', overwrite: 'true', pattern: '*.err'

  memory '5 GB'

  input:
    tuple val(sample_name), path(fq1), path(fq2)

  output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: checkValidity_fqs
    path("${sample_name}.err", emit: checkValidity_log)

  script:
    error_log = "${sample_name}.err"

    """
    is_ok=\$(fqtools validate ${fq1} ${fq2})

    if [ \$is_ok == 'OK' ]; then printf 'OK' && printf "" >> ${error_log}; else echo "error: sample did not pass fqtools validation check\n \$is_ok" >> ${error_log}; fi
    """

  stub:
    error_log  = "${sample_name}.err"

    """
    printf OK
    touch ${error_log}
    """
}


process countReads {
  /**
  * fail sample if there are < 1M raw reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*.err'

  memory '5 GB'

  input:
    tuple val(sample_name), path(fq1), path(fq2), val(is_ok)
    val(read_n_threshold)

  when:
    is_ok == 'OK'

  output:
    tuple val(sample_name), path(fq1), path(fq2), stdout, emit: countReads_fqs
    path("${sample_name}.err", emit: countReads_log)

  script:
    error_log = "${sample_name}.err"

    """
    num_reads=\$(fqtools count ${fq1} ${fq2})

    if (( \$num_reads > ${read_n_threshold} )); then printf "" >> ${error_log} && printf "${sample_name}"; else echo "error: sample did not have >= ${read_n_threshold} pairs of raw reads (it only contained \$num_reads)" >> ${error_log} && printf "fail"; fi
    """

  stub:
    error_log = "${sample_name}.err"

    """
    fqtools -v
    printf ${sample_name}
    touch ${error_log}
    """
}


process fastp {
  /**
  * confirm that there > 1M reads after cleaning with fastp
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/raw_read_QC_reports", mode: 'copy', pattern: '*.json'
  publishDir "${params.output_dir}/${sample_name}", mode: 'copy', overwrite: 'true', pattern: '*.err'

  memory '5 GB'

  input:
    tuple val(sample_name), path(fq1), path(fq2), val(run_fastp)
    val(read_n_threshold)

  when:
    run_fastp =~ /${sample_name}/

  output:
    tuple val(sample_name), path("${sample_name}_cleaned_1.fq.gz"), path("${sample_name}_cleaned_2.fq.gz"), stdout, emit: fastp_fqs
    path("${sample_name}_fastp.json", emit: fastp_json)
    path("${sample_name}.err", emit: fastp_log)

  script:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}.err"

    """
    fastp -i ${fq1} -I ${fq2} -o ${clean_fq1} -O ${clean_fq2} -j ${fastp_json} -h ${fastp_html} --length_required 50 --average_qual 10 --low_complexity_filter --correction --cut_right --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20

    rm -rf ${fastp_html}

    num_reads=\$(fqtools count ${clean_fq1} ${clean_fq2})

    if (( \$num_reads > ${read_n_threshold} )); then printf "" >> ${error_log} && printf "${sample_name}"; else echo "error: after fastp, sample did not have >= ${read_n_threshold} pairs of reads (it only contained \$num_reads)" >> ${error_log} && printf "fail"; fi
    """

  stub:
    clean_fq1  = "${sample_name}_cleaned_1.fq.gz"
    clean_fq2  = "${sample_name}_cleaned_2.fq.gz"
    fastp_json = "${sample_name}_fastp.json"
    fastp_html = "${sample_name}_fastp.html"
    error_log  = "${sample_name}.err"

    """
    fastp --version
    fqtools -v
    printf ${sample_name}
    touch ${error_log}
    touch ${clean_fq1}
    touch ${clean_fq2}
    touch ${fastp_json}
    touch ${fastp_html}
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
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)

  output:
    path("*", emit: fastQC_report)

  when:
    enough_reads =~ /${sample_name}/

  script:
    fqc_log = "${sample_name}fqc.log"

    """
    cat ${fq1} ${fq2} > ${sample_name}.fq.gz
    fastqc ${sample_name}.fq.gz 1> ${fqc_log}
    rm ${sample_name}.fq.gz
    """

  stub:
    report="${sample_name}.html"

    """
    fastqc -v
    touch ${report}
    """
}


process multiQC {
  /**
  * run MultiQC to aggregate fastQC reports
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
    multiqc --version
    touch ${multiQC_report}
    """
}


process trimGalore {
  /**
  * run trimGalore on fastq pair
  */
  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/preprocessing/Trim", mode: 'copy', pattern: '*.trimming_report.txt'

  memory '10 GB'

  input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)

  output:
    tuple path("${fq1}_trimming_report.txt"), path("${fq2}_trimming_report.txt"), emit: tg_reports
    tuple val(sample_name), path("*_val_1.*"), path("*_val_2.*"), emit: tg_fqs

  script:
    tg_log = "${sample_name}_tg.log"
    
    """
    trim_galore --paired ${fq1} ${fq2} 1> ${tg_log}
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

  cpus 4
  memory '15 GB'

  publishDir "${params.output_dir}/${sample_name}/preprocessing/Map", mode: 'copy', pattern: '*_alnStats.txt'
  publishDir "${params.output_dir}/phylo", mode: 'copy', pattern: '*.missing.bed'

  input:
    tuple val(sample_name), path(fq1), path(fq2), val(enough_reads)
    path(ref_bt2index)

  output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: bam
    path("${sample_name}.missing.bed"), emit: missing_bam
    path("${sample_name}_alnStats.txt"), emit: map_stats

  script:
    bam = "${sample_name}.sorted.bam"
    bai = "${sample_name}.bam.bai"
    stats = "${sample_name}_alnStats.txt"
    scripts = "${baseDir}/bin"

    """
    echo ${ref_bt2index}
    bowtie2 --very-sensitive -p ${task.cpus} -x ${ref_bt2index}/${params.ref} -1 ${fq1} -2 ${fq2} 2> ${sample_name}_alnStats.txt | samtools view -h - | samtools sort - -o ${bam}
    samtools depth ${bam} -aa | python3 ${scripts}/filter_depth.py 5 > ${sample_name}.missing.bed
    """

  stub:
    bam = "${sample_name}.sorted.bam"
    missing_bed = "${sample_name}.missing.bed"
    stats_txt = "${sample_name}_alnStats.txt"

    """
    python3 -V
    samtools --version
    bowtie2 --version
    touch ${bam}
    touch ${missing_bed}
    touch ${stats_txt}
    """
}


process picard {
  /**
  * run picard to perform deduplication on reads
  */
  tag { sample_name }

  memory '5 GB'

  publishDir "${params.output_dir}/bams", mode: 'copy', pattern: '*bam'

  input:
    tuple val(sample_name), path(bam)

  output:
    tuple val(sample_name), path("${sample_name}_grouped.bam"), emit: grouped_bam

  script:
    dedup_log = "${sample_name}_dedup.log"
    dedup_line="  java -jar /usr/local/bin/picard.jar MarkDuplicates INPUT=${bam} OUTPUT=${sample_name}_dedup.bam METRICS_FILE=${sample_name}_metrics.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true 1> ${dedup_log}"
    
    """
    java -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups I=${bam} O=${sample_name}_grouped.bam SORT_ORDER=coordinate RGID=1 RGPU=bc RGLB=lib RGPL=illumina RGSM=${sample_name} CREATE_INDEX=True
    """

  stub:
    grouped_bam = "${sample_name}_grouped.bam"
    metrics = "${sample_name}_metrics.txt"

    """
    if [ ! -f /usr/local/bin/picard.jar ]; then
      echo "ERROR: Picard JAR not found at /usr/local/bin/picard.jar for stub run." >&2
      exit 1
    fi
    echo "Picard JAR found at /usr/local/bin/picard.jar (stub check)."
    
    touch ${grouped_bam}
    touch ${metrics}
    """
}


process summarise {
  /**
  * extract mapping stats and produce JSON summmary
  */
  tag { sample_name }

  publishDir "${params.output_dir}/mapping_stats", mode: 'copy', pattern: '*.json'

  memory '5 GB'

  input:
    tuple val(sample_name), path(bam)

  output:
    path "${sample_name}.GG", emit: GG
    tuple val(sample_name), path(bam), path("${sample_name}_mapstats.json"), emit: bam_pre_out

  script:
    """
    samtools depth -a ${bam} > ${sample_name}.doc.bed
    python3 ${baseDir}/bin/gini.py ${sample_name}.doc.bed -G 5 1000 50 > ${sample_name}.GG
    samtools index ${bam}
    samtools stats ${bam} > tmp.stats
    python3 ${baseDir}/bin/parse_samtools_stats.py ${bam} tmp.stats ${sample_name}.doc.bed ${sample_name}.GG > ${sample_name}_mapstats.json
    """

  stub:
    gg_file = "${sample_name}.GG"
    doc_bed = "${sample_name}.doc.bed"
    mapstats_json = "${sample_name}_mapstats.json"
    """
    samtools --version
    python3 -V
    touch ${gg_file}
    touch ${doc_bed}
    touch ${mapstats_json}
    """
}
