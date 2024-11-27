// modules for downloading and preprocessing data


process getRefData {
  /**
  * download genome and annotation files
  */

  publishDir "${params.output_dir}/REFDATA/", mode: 'copy', overwrite: 'true', pattern: '*{.fasta,.gff}'

  memory '5 GB'

  input:
    val(genome_id)

  output:
    tuple path("${genome_id}.fasta"), path("${genome_id}.gff"), emit: refdata


  script:
    """
    cp ${baseDir}/resources/ref/${genome_id}/* ./
    """

  stub:
    fasta = "${genome_id}.fasta"
    gff = "${genome_id}.gff"

    """
    echo ">chr0
    ATCG
    >chr1
    ATCG
    >chr2
    ATCG
    >chr3
    ATCG
    >chr4
    ATCG" >> ${fasta}
    touch $gff
    """
}


process indexRefData {
  /**
  * index and preprocess reference fasta files
  */

  memory '5 GB'

  input:
    tuple path(fasta), path(gff)
    val(genome_id)

  output:
    path("./${genome_id}"), emit: bt2_index
    path("${fasta}.fai"), emit: reffaidx

  script:
    """
    samtools faidx ${fasta}
    mkdir ./${genome_id} && cd ./${genome_id}
    bowtie2-build ../${fasta} ${genome_id} && cd ../
    """

  stub:
    bt2_index = "./${genome_id}"
    reffaidx = "${fasta}.fai"

    """
    mkdir ${bt2_index}
    touch ${reffaidx}
    touch ${genome_id}
    """
}
