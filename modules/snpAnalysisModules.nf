// modules for identifying SNPs in mapped reads

process getMOI {
  /**
  * Investigate sample heterozygosity
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
  publishDir "${params.output_dir}/${sample_name}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.json'
  publishDir "${params.output_dir}/${sample_name}/MOI", mode: 'copy', overwrite: 'true', pattern: '*.png'

  memory '5 GB'

  input:
  tuple val(sample_name), path(grouped_bam)
  tuple path(fasta), path(gff)
  val(ref_id)

  output:
  path("${sample_name}.heterozygosity.json"), emit: moi_json
  tuple path("${sample_name}.bafscore.png"), path("${sample_name}.bafplot.png"), emit: moi_plots


  script:
  scripts = "${workflow.launchDir}/bin"
  """
  samtools faidx ${fasta}
  gatk --java-options "-Xmx8G" CreateSequenceDictionary R=${fasta} O=${ref_id}.dict
  gatk --java-options "-Xmx8G" HaplotypeCaller -R ${fasta} -I ${grouped_bam} -O ${sample_name}.vcf
  python3 ${scripts}/get_moi.py ${sample_name}.vcf ${sample_name}
  """

  stub:
  """
  touch ${sample_name}.heterozygosity.json
  touch ${sample_name}.bafscore.png
  touch ${sample_name}.bafplot.png
  touch ${sample_name}.vcf
  """
}

process phylo {
  /**
  * Run Phylogenetic analysis for sequence and SNP typing schemes
  */

  publishDir "${params.output_dir}/${sample_name}/phylo", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  path(grouped_bam)
  tuple path(fasta), path(gff)
  val(database)
  path(vars_yaml)

  output:
  path("*.gvcf.gz"), emit: st_vcf
  path("snps_merged.vcf.gz"), emit: snp_vcf
  path("snp_tree.nwk"), emit: snp_nwk
  path("*png"), emit: snp_png
  path("iqtree.json"), emit: iqtree_json
  path("*iqtree"), emit: iqtree
  path("*mldist"), emit: mldist

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  python3 ${scripts}/get_MLVA.py --bams ./ ${vars_yaml} ${fasta}
  a=`ls -la *vcf | wc -l`
  if [ \$a -gt 2 ]
    then
      bcftools merge --force-samples -g ${fasta} -m snps ./*gvcf.gz -o snps.tmp.gvcf
      bcftools filter -O z -o snps_merged.vcf.gz -i '%QUAL>=30' snps.tmp.gvcf
      Rscript ${scripts}/snpPhylo.R -v snps_merged.vcf.gz
    else
      touch NO_TREE.png
  fi
  """

  stub:
  """
  """
}

process findSNPs {
  /**
  * Call snps in mapped reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/GVCFs", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'

  memory '8 GB'

  input:
  tuple val(sample_name), path(grouped_bam)
  tuple path(fasta), path(gff)
  val(ref_id)

  output:
  path("${sample_name}.vcf.gz"), emit: gvcf
  path("*MERGE_AND_PLOT"), emit: merge_and_plot

  script:
  error_log = "${sample_name}.err"

  """
  bcftools mpileup -Ov --gvcf 0 -f ${fasta} ${grouped_bam} | bcftools call --skip-variants indels -m --gvcf 0 -o ${sample_name}.gvcf
  bgzip -ci ${sample_name}.gvcf > ${sample_name}.vcf.gz

  a=`ls -la *vcf | wc -l`
  if [ \$a -gt 2 ]
    then
      touch ./MERGE_AND_PLOT;
    else
      touch ./DONT_MERGE_AND_PLOT;
  fi
  """

  stub:
  error_log  = "${sample_name}.err"

  """
  touch ${sample_name}.vcf.gz
  touch ${sample_name}.vcf
  touch ./MERGE_AND_PLOT
  """
}

process makeVCFtree {
  /**
  * Plot SNPs across each chromosome
  */

  publishDir "${params.output_dir}/GVCFs", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
  publishDir "${params.output_dir}/SNP_phylos", mode: 'copy', overwrite: 'true', pattern: '*.png'
  publishDir "${params.output_dir}/SNP_phylos", mode: 'copy', overwrite: 'true', pattern: '*.nwk'

  memory '5 GB'

  input:
  path(gvcf_files)
  val(database)
  tuple path(fasta), path(gff)
  path(merge_and_plot)

  output:
  path("merged.vcf.gz"), emit: vcf
  path("*.png"), emit: snp_phylo_pngs
  path("*.nwk"), emit: newick_tree

  when:
  merge_and_plot = /MERGE\_AND\_PLOT/

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  ln -s ${database}
  if [ -d ${database} ]; then
    cp ${database}/*.vcf.gz ./;
  fi

  for x in ./*vcf.gz; do
    bcftools index \$x;
  done

  a=`ls -la *vcf.gz | wc -l`
  if [ \$a -gt 1 ]
    then
      bcftools merge --force-samples -g ${fasta} -o merged.vcf ./*vcf.gz
      gzip merged.vcf
    else
      mv *vcf.gz merged.vcf.gz
      bcftools index merged.vcf.gz
  fi
  if [ \$a -gt 2 ]
    then
      Rscript ${scripts}/snpPhylo.R -v merged.vcf.gz
    else
    touch NO_TREE.png
    touch NO_TREE.nwk
  fi
  """

  stub:
  """
  ln -s ${database}
  if [ -d ${database} ]; then
    cp ${database}/*.vcf.gz ./;
  fi

  for x in ./*vcf.gz; do
    touch \${x}.png;
  done
  touch merged.vcf.gz
  touch tree.nwk
  """
}

process makeJSON {
  /**
  * Make JSON for SNP module
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(grouped_bam)
  path(mapstats_json)
  path(snp_nwk)
  path(iqtree_json)

  output:
  path("${sample_name}_SNPreport.json"), emit: snpreport_json

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  python3 ${scripts}/make_SNP_json.py ${mapstats_json} ${snp_nwk} ${iqtree_json} > ${sample_name}_SNPreport.json
  """

  stub:
  """
  touch ${sample_name}_SNPreport.json
  """
}

process makeReport {
  /**
  * Make report for SNP module
  */

  tag { sample_name }

  publishDir "${params.output_dir}/${sample_name}/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(snpreport_json)
  path(snpQC_plots)

  output:
  path("${sample_name}_report.pdf"), emit: report_pdf

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  python3 ${scripts}/makeReport.py ${snpreport_json} ${snpQC_plots} ${sample_name}
  """

  stub:
  """
  touch ${sample_name}_report.pdf
  """
}
