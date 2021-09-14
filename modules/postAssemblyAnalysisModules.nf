// modules for post assembly genome analysis

process makeChromosomeFastas {
  /**
  * generates fastas of individual chromosomes from a set of assembly fastas scaffolded using ABACAS
  */

  memory '5 GB'

  input:
  tuple val(assemblies)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("*fasta"), emit: chromosome_fastas

  script:
  """
  echo ${assemblies}
  python3 ${baseDir}/scripts/makeChromosomeFastas.py ${assemblies}
  """

}

process chromosomeAlignment {
  /**
  * Run spades assembly on preprocessed reads
  */
  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/Assembly", mode: 'copy', overwrite: 'true'

  cpus 8
  memory '15 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)

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
