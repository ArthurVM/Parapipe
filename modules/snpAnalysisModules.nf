// modules for identifying SNPs in mapped reads

process getMOI {
  /**
  * Investigate sample heterozygosity
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'
  publishDir "${params.output_dir}/${sample_name}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.json'
  publishDir "${params.output_dir}/${sample_name}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.png'

  memory '15 GB'
  cpus 2

  input:
  tuple val(sample_name), path(bam), path(mapstats)
  tuple path(fasta), path(gff)
  val(ref_id)

  output:
  tuple val(sample_name), path(bam), path(mapstats), path("${sample_name}.heterozygosity.json"), path("${sample_name}.*.png"), path("${sample_name}.vcf.gz"), emit: bam_pp_moi
  // path("${sample_name}.heterozygosity.json"), emit: moi_json
  // tuple path("${sample_name}.bafscore.png"), path("${sample_name}.bafplot.png"), emit: moi_plots

  script:
  scripts = "${baseDir}/bin"
  """
  freebayes -p 2 -P 0 -C 2 -F 0.05 --min-coverage 5 --min-repeat-entropy 1.0 --min-alternate-qsum 30 -q 13 -m 1 --strict-vcf -f ${fasta} ${bam} > ${sample_name}.tmp.vcf
  sed -i \'s/##FORMAT=<ID=GQ,Number=1,Type=Integer/##FORMAT=<ID=GQ,Number=1,Type=String/\' ${sample_name}.tmp.vcf
  sed -i \'s/unknown/${sample_name}/\' ${sample_name}.tmp.vcf
  python3 ${scripts}/filter_vcf.py  ${sample_name}.tmp.vcf --qual 30 --only_snp -o ${sample_name}.vcf
  python3 ${scripts}/get_moi.py ${sample_name}.vcf ${sample_name}

  bgzip -c ${sample_name}.vcf > ${sample_name}.vcf.gz
  """

  stub:
  """
  freebayes --version
  touch ${sample_name}.heterozygosity.json
  touch ${sample_name}.bafscore.png
  touch ${sample_name}.bafplot.png
  touch ${sample_name}.vcf.gz
  """
}

process phylo {
  /**
  * Run Phylogenetic analysis for sequence and SNP typing schemes
  */

  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.json'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.nxs'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.png'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.csv'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.nwk'
  publishDir "${params.output_dir}/phylo", mode: 'copy', overwrite: 'true', pattern: '*.snps.bed'

  memory '15 GB'
  cpus 8

  input:
  path(bam)
  path(mapstats_jsons)
  path(vcf)
  tuple path(fasta), path(gff)
  val(database)
  path(vars_yaml)

  output:
  path("snp_tree.nwk"), emit: snp_nwk
  path("*_ST_stats.csv"), emit: st_stats
  path("*.snps.bed"), emit: snps_bed
  path("allele_matrix.csv"), emit: allele_matrix
  path("*.png"), emit: snp_png
  path("iqtree.json"), emit: iqtree_json
  path("*.nxs"), emit: nxs_aln
  path("allele_stats.json"), emit: al_stats_json

  script:
  scripts = "${baseDir}/bin"
  """
  for v in ./*vcf.gz; do
      tabix -f -p vcf \$v
  done

  python3 ${scripts}/get_MLVA_consensus.py ${vars_yaml} --bams ./*bam

  a=`ls -la *vcf.gz | wc -l`
  if [ \$a -gt 2 ]
    then
      python3 ${scripts}/vcf_to_allelematrix.py --vcfs *vcf.gz --mapstats ${mapstats_jsons}
      python3 ${scripts}/make_WG_tree.py allele_matrix.csv
    else
      touch NO_SNP_TREE.png
  fi

  touch snp_tree.nwk
  """

  stub:
  """
  touch snp_tree.nwk
  touch stub_ST_stats.csv
  touch stub.snps.bed
  touch allele_matrix.csv
  touch stub.png
  touch iqtree.json
  touch stub.nxs
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
  tuple val(sample_name), path(bam), path(mapstats_json), path(moi_json), path(moi_png), path(vcf)
  path(snp_nwk)
  path(iqtree_json)
  path(st_stats)
  path(al_stats_json)

  output:
  path("${sample_name}_report.json"), emit: report_json

  script:
  scripts = "${baseDir}/bin"
  st_stats_csv = "${sample_name}_ST_stats.csv"
  """
  python3 ${scripts}/make_report_json.py ${mapstats_json} ${iqtree_json} ${st_stats_csv} ${al_stats_json} > ${sample_name}_report.json
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
  tuple val(sample_name), path(bam), path(mapstats_json), path(moi_json), path(moi_png), path(vcf)
  path(env_json)
  path(snpreport_json)
  path(snp_png)

  output:
  path("${sample_name}_report.pdf"), emit: report_pdf

  script:
  scripts = "${baseDir}/bin"
  """
  python3 ${scripts}/make_sample_report_pdf.py --env ${env_json} --report ${snpreport_json} --id ${sample_name} --moi_json ${moi_json} --moi_plots ${moi_png} --png ${snp_png}
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

  output:
  path("Parapipe_report.pdf"), emit: run_report_pdf

  script:
  scripts = "${baseDir}/bin"
  """
  python3 ${scripts}/make_run_report_pdf.py --env ${env_json} --report ${snpreport_jsons} --moi_json ${moi_json}
  """

  stub:
  """
  touch Parapipe_report.pdf
  """
}
