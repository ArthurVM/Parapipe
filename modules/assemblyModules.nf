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

  publishDir "${params.output_dir}/$sample_name/Assembly", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(fasta), path(gff), path(cds), path(gaf)
  path(spades_assembly)

  output:
  tuple path("aligned_stats"), path("basic_stats"), path("genome_stats")
  path("contigs_reports")
  tuple path("icarus.html"), path("report.html")
  path("quast.log")
  tuple path("*.tex")
  tuple path("*.tsv")
  tuple path("*.txt")
  path("quast.log")

  script:
  """
  quast.py -o ./ -r ${fasta} -g ${gff } ${spades_assembly}
  """

  stub:
  """
  """
}
