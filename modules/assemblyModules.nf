// modules for assembly and annotation

process spades {
  /**
  * Run spades assembly on preprocessed reads
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/Assembly", mode: 'copy', overwrite: 'true'

  cpus 8
  memory '15 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(trimmed_fq1), path(trimmed_fq2)

  output:
  path("scaffolds.fasta"), emit: scaffolds
  path("spades.log"), emit: log

  script:
  """
  spades.py -1 ${trimmed_fq1} -2 ${trimmed_fq2} --careful --cov-cutoff auto -k auto --threads ${task.cpus} --memory 15 -o ./
  """

  stub:
  """
  touch scaffolds.fasta
  touch spades.log
  """
}

process quast {
  /**
  * Run quast on spades assemblies
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/Assembly/quast/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(fasta), path(gff), path(cds), path(gaf)
  path(spades_assembly)

  output:
  path("icarus.html"), emit: icarus
  path("quast.log"), emit: log
  tuple path("report.html"), path("report.tex"), path("report.tsv"), path("report.txt"), emit: reports
  tuple path("transposed_report.tex"), path("transposed_report.tsv"), path("transposed_report.txt"), emit: transposed_reports
  tuple path("aligned_stats"), path("basic_stats"), path("genome_stats")
  path("contigs_reports")

  script:
  """
  quast.py -o ./ -r ${fasta} -g ${gff } ${spades_assembly}
  """

  stub:
  """
  touch icarus.html
  touch quast.log
  touch report.{html,tex,tsv,txt}
  touch transposed_report.{tex,tsv,txt}
  mkdir {aligned,basic,genome}_stats
  """
}

process indexAssembly {
  /**
  * index a genome assembly
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/Assembly/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(scaffolds)

  output:
  path("${sample_name}"), emit: bt2_index

  script:
  """
  mkdir ./${sample_name} && cd ./${sample_name}
  bowtie2-build ../${scaffolds} ${sample_name} && cd ../
  """

  stub:
  bt2_index = "./${params.genome}"

  """
  mkdir ${bt2_index}
  """
}

process map2SPAdesFasta {
  /**
  * map reads to their genome assembly using Bowtie2
  */

  tag { sample_name }

  cpus 8
  memory '15 GB'

  publishDir "${params.output_dir}/$sample_name/Assembly/Map", mode: 'copy'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(trimmed_fq1), path(trimmed_fq2)
  path(ref_bt2index)

  output:
  path("${sample_name}.sorted.bam"), emit: bam
  path("${sample_name}_alnStats.txt"), emit: map_stats


  script:
  bam = "${sample_name}.sorted.bam"
  stats = "${sample_name}_alnStats.txt"

  """
  echo $ref_bt2index
  bowtie2 -p ${task.cpus} -x ${ref_bt2index}/${ref_bt2index} -1 $fq1 -2 $fq2 2> ${sample_name}_alnStats.txt | samtools view -h - | samtools sort - -o ${bam}
  samtools index ${bam}
  """

  stub:
  bam = "${sample_name}.sorted.bam"
  stats_txt = "${sample_name}_alnStats.txt"

  """
  touch ${bam}
  touch ${stats_txt}
  """
}

process pilon {
  /**
  * map reads to their genome assembly using Bowtie2
  * TODO: this throws weird errors, needs fixing
  */

  tag { sample_name }

  errorStrategy = {task.attempt <= 3 ? 'retry' : 'ignore'}
  memory = {5 GB * task.attempt}
  cpus 4
  maxRetries = 3

  publishDir "${params.output_dir}/$sample_name/Assembly/Pilon", mode: 'copy'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(bam)
  path(scaffolds)

  output:
  path("${sample_name}.fasta"), emit: pilon_fasta
  path("${sample_name}.vcf"), emit: pilon_vcf
  path("pilon.log")
  path("pilon.err")

  script:
  pilon_stdout = "pilon.log"
  pilon_stderr = "pilon.err"
  """
  java -jar /usr/local/bin/pilon-1.24.jar --genome ${scaffolds} --bam ${bam} --output ${sample_name} --vcf > ${pilon_stdout} 2> ${pilon_stderr}
  """

  stub:
  """
  touch ${sample_name}.fasta
  touch ${sample_name}.vcf
  """
}

process abacas {
  /**
  * Map scaffolds to reference using ABACAS
  */

  tag { sample_name }

  memory = '10 GB'
  cpus 4

  publishDir "${params.output_dir}/$sample_name/Assembly/ABACAS", mode: 'copy'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(scaffolds)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("${sample_name}.abacas.fasta"), emit: abacas_fasta
  path("${sample_name}*.features.tab"), emit: features
  path("unused_contigs.out"), emit: unused_contigs

  script:
  mummer_module = "nucmer"
  ref_union="${fasta}.union"
  abacas_fasta="${scaffolds}_${fasta}.union.fasta"
  crunch="${scaffolds}_${fasta}.union.fasta.crunch"
  tab="${scaffolds}_${fasta}.union.fasta.tab"
  """
  perl \$ABACAS_DIR/joinMultifasta.pl ${fasta} ${fasta}.union
  perl \$ABACAS_DIR/abacas.1.3.1.pl -r ${ref_union} -q ${scaffolds} -p ${mummer_module}
  perl \$ABACAS_DIR/splitABACASunion.pl ${fasta} ${ref_union} ${abacas_fasta} ${crunch} ${tab} ${sample_name}
  """

  stub:
  """
  """
}
