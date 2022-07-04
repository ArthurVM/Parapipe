// modules for identifying SNPs in mapped reads

process findSNPs {
  /**
  * Call snps in mapped reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true', pattern: '*.var.vcf'
  publishDir "${params.output_dir}/GVCFs", mode: 'copy', overwrite: 'true', pattern: '*.vcf.gz'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(dedup_bam)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("${sample_name}.var.vcf"), emit: vcf
  path("${sample_name}.vcf.gz"), emit: gvcf
  path("*MERGE_AND_PLOT"), emit: merge_and_plot

  script:
  error_log = "${sample_name}.err"

  """
  bcftools mpileup -Ov --gvcf 0 -f ${fasta} ${dedup_bam} | bcftools call -m --gvcf 0 -o ${sample_name}.gvcf
  bgzip -ci ${sample_name}.gvcf > ${sample_name}.vcf.gz
  samtools mpileup --skip-indels --BCF -f ${fasta} ${dedup_bam} | bcftools call --skip-variants indels -m -O v --variants-only -o ${sample_name}.var.vcf -

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
  tuple path(fasta), path(gff), path(cds), path(gaf)
  path(merge_and_plot)

  output:
  path("merged.vcf"), emit: vcf
  path("*.png"), emit: snp_phylo_pngs
  path("*.nwk"), emit: newick_tree

  when:
  merge_and_plot = /MERGE\_AND\_PLOT/

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  for x in ./*vcf.gz; do
    bcftools index \$x;
  done
    bcftools merge -g ${fasta} -o merged.vcf ./*vcf.gz
  Rscript ${scripts}/snpPhylo.R -v merged.vcf
  """

  stub:
  """
  touch merged.vcf
  for x in ./*vcf.gz; do
    touch \${x}.png;
  done
  """
}

process preprocessForPlotting {
  /**
  * preprocesses FASTA and VCF files for plotting
  */

  tag { sample_name }

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("tmp.descStripped.fasta"), emit: preprocessed_fasta
  path("vcf_split"), emit: vcf_path

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  echo "Preprocessing ${fasta}..."
  python3 ${scripts}/removeFASTAdesc.py ${fasta}
  echo "Splitting ${vcf} by chromosome ID..."
  python3 ${scripts}/splitVCF.py ${vcf}
  """

  stub:
  """
  touch tmp.descStripped.fasta
  mkdir vcf_split
  """
}

process plotSNP {
  /**
  * Plot SNPs across each chromosome
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(vcf_path)
  path(preprocessed_fasta)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("SNPDist.pdf"), emit: pdf

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  echo "Creating plots..."
  Rscript ${scripts}/SNPdist.R -v ${vcf_path} -f ${preprocessed_fasta} -g ${gff}
  """

  stub:
  """
  touch SNPDist.pdf
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
  path(mapstats_json)
  path(newick_tree)

  output:
  path("${sample_name}_SNPreport.json"), emit: snpreport_json

  script:
  scripts = "${workflow.launchDir}/scripts"
  """
  python3 ${scripts}/makeReport.py ${mapstats_json} ${newick_tree} > ${sample_name}_SNPreport.json
  """

  stub:
  """
  touch SNPDist.pdf
  """
}
