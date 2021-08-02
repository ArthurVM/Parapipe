// modules for downloading and preprocessing data

process downloadRefData {
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

  stub:
  fasta = "${genome_id}.fasta"
  gff = "${genome_id}.gff"
  cds_fasta = "${genome_id}_cds.fasta"

  """
  touch $fasta
  touch $gff
  touch $cds_fasta
  """

}

process indexRefData {
  /**
  * index and preprocess reference fasta files
  */

  publishDir "${params.output_dir}/REFDATA/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple path(fasta), path(gff), path(cds)
  val(genome_id)

  output:
  path("./${params.genome}"), emit: bt2_index
  path("${fasta}.fai"), emit: reffaidx
  path("${genome_id}.dict"), emit: refdict

  script:
  """
  samtools faidx ${fasta}
  mkdir ./${genome_id} && cd ./${genome_id}
  bowtie2-build ../${fasta} ${genome_id} && cd ../
  java -jar /usr/local/bin/picard.jar CreateSequenceDictionary REFERENCE=${fasta} OUTPUT=${genome_id}.dict
  """

  stub:
  bt2_index = "./${params.genome}"
  reffaidx = "${fasta}.fai"
  refdict = "${genome_id}.dict"

  """
  mkdir ${bt2_index}
  touch ${reffadx}
  touch ${refdict}
  """
}
