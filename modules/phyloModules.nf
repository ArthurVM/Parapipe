// modules for identifying SNPs in mapped reads

process callVariants {
  /**
  * Call variants
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/VCF", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'

  memory '15 GB'
  cpus 2

  input:
  tuple val(sample_name), path(bam), path(mapstats)
  tuple path(fasta), path(gff)
  val(ref_id)

  output:
  tuple val(sample_name), path(bam), path(mapstats), path("${sample_name}.vcf.gz"), emit: bam_vcf_tup
  path("${sample_name}.vcf.gz"), emit: vcf


  script:
  scripts = "${baseDir}/bin"
  """
  freebayes -p 2 -P 0 -C 2 -F 0.05 --min-coverage 5 --min-repeat-entropy 1.0 --min-alternate-qsum 30 -q 13 -m 1 --strict-vcf -f ${fasta} ${bam} > ${sample_name}.tmp.vcf
  sed -i \'s/##FORMAT=<ID=GQ,Number=1,Type=Integer/##FORMAT=<ID=GQ,Number=1,Type=String/\' ${sample_name}.tmp.vcf
  sed -i \'s/unknown/${sample_name}/\' ${sample_name}.tmp.vcf

  python3 ${scripts}/filter_vcf.py  ${sample_name}.tmp.vcf --qual 30 --only_snp -o ${sample_name}.vcf

  bgzip -c ${sample_name}.vcf > ${sample_name}.vcf.gz
  """

  stub:
  """
  touch ${sample_name}.vcf.gz
  """
}

process moi {
  /**
  * Investigate multiplicity of infection
  */

  publishDir "${params.output_dir}/MOI", mode: 'copy', overwrite: 'true', pattern: 'filtered_merged.vcf.gz'
  publishDir "${params.output_dir}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.json'
  publishDir "${params.output_dir}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.png'

  memory '15 GB'
  cpus 2

  input:
  path(vcfs)
  tuple path(fasta), path(gff)

  output:
  path("*.moi.json"), emit: moi_json
  path("*.png"), emit: moi_png
  path("fws.csv"), emit: fws_csv
  path("filtered_merged.vcf.gz"), emit: merged_vcf

  script:
  scripts = "${baseDir}/bin"
  """
  for v in ./*vcf.gz; do
      bcftools index \$v
  done

  bcftools merge -g ${fasta} -m snps ./*vcf.gz > merged.vcf

  vcftools --vcf ./merged.vcf --remove-indels --maf 0.05 --max-missing 0.1 --minQ 30 --min-meanDP 5 --minDP 5 --recode --stdout > filtered_merged.vcf

  Rscript ${scripts}/getFWS.R --vcf=filtered_merged.vcf

  python3 ${scripts}/get_moi.py filtered_merged.vcf fws.csv
  
  bgzip filtered_merged.vcf
  """

  stub:
  """
  touch fws.csv
  touch filtered_merged.vcf.gz
  touch test.moi.json
  touch test.png
  """
}

process molTyping_phylo {
  /**
  * Run Phylogenetic analysis for molecular typing schemes
  */

  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: 'typing_report.json'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.nxs'

  memory '15 GB'
  cpus 8

  input:
  path(bam)
  path(vars_yaml)

  output:
  path("typing_report.json"), emit: typingreport_json
  path("*.nxs"), emit: nxs_aln

  script:
  scripts = "${baseDir}/bin"
  """
  python3 ${scripts}/get_MLVA_consensus.py ${vars_yaml} --bams ./*bam
  """

  stub:
  """
  touch iqtree.json
  touch stub.nxs
  touch allele_stats.json
  """
}

process wgSNV_phylo {
  /**
  * Run Phylogenetic analysis for SNP typing scheme
  */

  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.json'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.png'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.snps.bed'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: 'allele_matrix.csv'

  memory '15 GB'
  cpus 8

  input:
  path(mapstats_jsons)
  path(vcf)
  tuple path(fasta), path(gff)
  path(database)
  val(mincov)

  output:
  path("*.snps.bed"), emit: snps_bed
  path("allele_matrix.csv"), emit: allele_matrix
  path("absolute_distance.csv"), emit: dist_matrix
  path("*.png"), emit: snp_png
  path("allele_stats.json"), emit: al_stats_json

  script:
  scripts = "${baseDir}/bin"
  """
  for v in ./*vcf.gz; do
      tabix -f -p vcf \$v
  done

  python3 ${scripts}/vcf_to_allelematrix.py --vcfs *vcf.gz --mapstats ${mapstats_jsons} --database ${database}
  python3 ${scripts}/make_WG_tree.py allele_matrix.csv ${mincov}
  """

  stub:
  """
  touch stub.snps.bed
  touch allele_matrix.csv
  touch stub.png
  touch allele_stats.json
  """
}

process makeJSON {
  /**
  * Make JSON for this module
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(bam), path(mapstats_json), path(vcf)
  path(typingreport_json)
  path(al_stats_json)

  output:
  path("${sample_name}_report.json"), emit: report_json

  script:
  scripts = "${baseDir}/bin"
  st_stats_csv = "${sample_name}_ST_stats.csv"
  """
  python3 ${scripts}/make_report_json.py ${sample_name} ${mapstats_json} ${typingreport_json} ${al_stats_json} > ${sample_name}_report.json
  """

  stub:
  """
  touch ${sample_name}_report.json
  """
}

process makeSampleReports {
  /**
  * Make PDF reports for each sample
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(bam), path(mapstats_json), path(vcf)
  path(env_json)
  path(snpreport_json)
  path(snp_png)
  path(moi_json)
  path(moi_png)

  output:
  path("${sample_name}_report.pdf"), emit: report_pdf

  script:
  scripts = "${baseDir}/bin"
  """
  python3 ${scripts}/make_sample_report_pdf.py --env ${env_json} --report ${snpreport_json} --id ${sample_name} --moi_json ${sample_name}.moi.json --png ${snp_png}
  """

  stub:
  """
  touch ${sample_name}_report.pdf
  """
}

process makeRunReport {
  /**
  * Make PDF report for Parapipe run
  */

  publishDir "${params.output_dir}/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  path(env_json)
  path(snpreport_jsons)
  path(snp_png)
  path(moi_json)
  path(moi_png)
  path(multiqc_report)
  path(dist_matrix)

  output:
  path("Parapipe_report.html"), emit: run_report_html
  path("run_results.csv"), emit: results_csv

  script:
  scripts = "${baseDir}/bin"
  """
  python3 ${scripts}/html_report.py --env ${env_json} --moi_json ${moi_json} --report ${snpreport_jsons} --dist_matrix ${dist_matrix} --multiqc ${baseDir}/${params.output_dir}/multiqc_report.html
  """

  stub:
  """
  touch Parapipe_report.pdf
  """
}
