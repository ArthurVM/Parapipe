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


process indexGp60DB {
  /**
  * index the gp60 database
  */

  memory '5 GB'

  output:
    path("./gp60_db"), emit: gp60_bt2_db

  script:
    gp60_db = "${baseDir}/resources/gp60_probes/gp60_db.fa"
    """
    bowtie2-build ${gp60_db} gp60_db
    mkdir ./gp60_db && mv *bt2 ./gp60_db
    """

  stub:
    """
    bowtie2 --version
    mkdir gp60_db
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
    """

  stub:
    reffaidx = "${fasta}.fai"

    """
    samtools --version
    touch ${reffaidx}
    touch ${genome_id}
    """
}
