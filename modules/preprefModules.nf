// modules for downloading and preprocessing data

process DownloadRefData {
  /**
  * download genome and annotation files
  */

  publishDir "${params.output_dir}/REFDATA/", mode: 'copy', overwrite: 'true', pattern: '*{.fasta,.gff}'

  memory '5 GB'

  input:
  val(genome_id)

  output:
  tuple path("${genome_id}.fasta"), path("${genome_id}.gff"), path("${genome_id}_cds.fasta"), emit: refdata

  script:
  """
  python3 ${baseDir}/scripts/DLandIndex.py $genome_id
  """

}

process IndexRefData {
  /**
  * index and preprocess reference fasta files
  */

  publishDir "${params.output_dir}/REFDATA/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple path(fasta), path(gff), path(cds)
  val(genome_id)

  output:
  path("${genome_id}"), emit: indexBT2
  path("${fasta}.fai"), emit: reffaidx
  path("${genome_id}.dict"), emit: refdict

  script:
  """
  samtools faidx ${fasta}
  mkdir ./${genome_id} && cd ./${genome_id}
  bowtie2-build ../${fasta} ${genome_id} && cd ../
  java -jar /usr/local/bin/picard.jar CreateSequenceDictionary REFERENCE=${fasta} OUTPUT=${genome_id}.dict
  """
}
