// modules for read mapping and variant calling via Snippy


process snippyCall {
  /**
  * Map reads and call variants with Snippy
  */

  tag { sample_name }

  cpus 4
  memory '15 GB'

  input:
    tuple val(sample_name), path(fq1), path(fq2)
    tuple path(fasta), path(gff)
    val(ref_id)

  output:
    tuple val(sample_name), path("${sample_name}/${sample_name}.bam"), path("${sample_name}/${sample_name}.vcf.gz"), emit: bam_vcf
    path("${sample_name}/${sample_name}.vcf.gz"), emit: vcf
    path("${sample_name}/${sample_name}.vcf.gz.tbi"), emit: vcf_index
    path("${sample_name}/${sample_name}.raw.vcf.gz"), emit: raw_vcf
    path("${sample_name}/${sample_name}.aligned.fa"), emit: aligned_fa
    path("${sample_name}/${sample_name}.consensus.fa"), emit: consensus_fa
    path("${sample_name}/${sample_name}.tab"), emit: snp_tab

  script:
    """
    snippy \\
      --outdir ${sample_name} \\
      --prefix ${sample_name} \\
      --ref ${fasta} \\
      --R1 ${fq1} \\
      --R2 ${fq2} \\
      --cpus ${task.cpus}

    bgzip -f ${sample_name}/${sample_name}.vcf
    bgzip -f ${sample_name}/${sample_name}.raw.vcf
    tabix -f -p vcf ${sample_name}/${sample_name}.vcf.gz
    """

  stub:
    """
    snippy --version
    bgzip --version
    tabix --version
    mkdir -p ${sample_name}
    touch ${sample_name}/${sample_name}.bam
    touch ${sample_name}/${sample_name}.vcf.gz
    touch ${sample_name}/${sample_name}.vcf.gz.tbi
    touch ${sample_name}/${sample_name}.raw.vcf
    touch ${sample_name}/${sample_name}.aligned.fa
    touch ${sample_name}/${sample_name}.consensus.fa
    touch ${sample_name}/${sample_name}.tab
    """
}
