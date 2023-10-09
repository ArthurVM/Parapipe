// modules for identifying SNPs in mapped reads

process getMOI {
  /**
  * Investigate sample heterozygosity
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
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

process findSNPs_OLD {
  /**
  * Call snps in mapped reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true', pattern: '*.vcf'
  publishDir "${params.output_dir}/GVCFs", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(dedup_bam)
  tuple path(fasta), path(gff)

  output:
  path("${sample_name}.var.vcf"), emit: vcf
  path("${sample_name}.vcf.gz"), emit: gvcf
  path("*MERGE_AND_PLOT"), emit: merge_and_plot

  script:
  error_log = "${sample_name}.err"

  """
  bcftools mpileup -Ov --gvcf 0 -f ${fasta} ${dedup_bam} | bcftools call --skip-variants indels -m --gvcf 0 -o ${sample_name}.gvcf
  bgzip -ci ${sample_name}.gvcf > ${sample_name}.vcf.gz
  bcftools mpileup --skip-indels -f ${fasta} ${dedup_bam} | bcftools call --skip-variants indels -m -O v --variants-only -o ${sample_name}.var.vcf -

  a=`ll *vcf.gz | wc -l`
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
  touch ${sample_name}.var.vcf
  touch ${sample_name}.vcf.gz
  touch ./MERGE_AND_PLOT
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
  tuple val(sample_name), path(grouped_bam), path("${sample_name}.vcf"), emit: vcf
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

process preprocessForPlotting {
  /**
  * preprocesses FASTA and VCF files for plotting
  */

  tag { sample_name }

  memory '5 GB'

  errorStrategy 'ignore'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf)
  tuple path(fasta), path(gff)

  output:
  path("fasta_split"), emit: fasta_path
  path("vcf_split"), emit: vcf_path

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  echo "Splitting ${fasta}..."
  python3 ${scripts}/preprocess_FASTA_for_plotting.py ${fasta}
  mkdir fasta_split
  mv tmp.*.fasta fasta_split

  echo "Splitting ${vcf} by chromosome ID..."
  python3 ${scripts}/splitVCF.py ${vcf}
  """

  stub:
  """
  mkdir fasta_split
  mkdir vcf_split
  """
}

process plotSNP {
  /**
  * Plot SNPs across each chromosome
  */

  tag { sample_name }

  memory '5 GB'

  errorStrategy 'ignore'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf_path)
  path(fasta_path)
  tuple path(fasta), path(gff)

  output:
  path("snpQC_plots"), emit: snpQC_plots

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  echo "Creating plots..."
  Rscript ${scripts}/SNPdist.R -v ${vcf_path} -f ${fasta_path} -g ${gff}
  mkdir snpQC_plots
  mv *jpeg snpQC_plots
  """

  stub:
  """
  touch test.jpeg
  mkdir snpQC_plots
  mv *jpeg snpQC_plot
  """
}

process makeJSON {
  /**
  * Make JSON for SNP module
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(grouped_bam)
  path(mapstats_json)
  path(newick_tree)

  output:
  path("${sample_name}_SNPreport.json"), emit: snpreport_json

  script:
  scripts = "${workflow.launchDir}/bin"
  """
  python3 ${scripts}/make_SNP_json.py ${mapstats_json} ${newick_tree} > ${sample_name}_SNPreport.json
  """

  stub:
  """
  touch ${sample_name}_SNPreport.json
  """
}

process makeSNPreport {
  /**
  * Make report for SNP module
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/", mode: 'copy', overwrite: 'true'

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
