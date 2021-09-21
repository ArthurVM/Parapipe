// modules for identifying SNPs in mapped reads

process findSNPs {
  /**
  * Call snps in mapped reads
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/SNP", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(dedup_bam)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("${sample_name}.var.vcf"), emit: vcf

  script:
  error_log = "${sample_name}.err"

  """
  samtools mpileup --skip-indels --BCF -f ${fasta} ${dedup_bam} | bcftools call --skip-variants indels -m -O v --variants-only -o ${sample_name}.var.vcf -
  """

  stub:
  error_log  = "${sample_name}.err"

  """
  touch ${sample_name}.var.vcf
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

process findSTRs {
  /**
  * Find Short Tandem Repeats using TRF
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/STR", mode: 'copy', overwrite: 'true'

  memory '12 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(scaffolds)

  output:
  tuple path("trf.stdout"), path("trf.stderr"), emit: trf_log
  path("${scaffolds}.2.7.7.80.10.50.500.dat"), emit: trf_dat

  script:
  trf_stdout = "trf.stdout"
  trf_stderr = "trf.stderr"
  // There is some weirdness going on with TRF, where it is successfully running but producing non 0 exit status
  // trap is used here as a temporary fix
  """
  trap 'if [[ \$? != 0 ]]; then echo OVERRIDE EXIT; exit 0; fi' EXIT
  trf4.10.0-rc.2.linux64.exe ${scaffolds} 2 7 7 80 10 50 500 -d -h 1> ${trf_stdout} 2> ${trf_stderr}
  """

  stub:
  trf_dat = "${scaffolds}.2.7.7.80.10.50.500.dat"
  """
  touch ${trf_dat}
  touch trf.{stdout,stderr}
  """
}

process indexScaffolds {
  /**
  * index the assembled scaffolds
  */

  tag { sample_name }

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(scaffolds)

  output:
  path("./${sample_name}/"), emit: scaffold_bt2index

  script:
  """
  mkdir ./${sample_name} && cd ./${sample_name}
  bowtie2-build ../${scaffolds} ${sample_name} && cd ../
  """

  stub:
  bt2_index = "./${sample_name}"

  """
  mkdir ${bt2_index}
  """
}

process map2Scaffolds {
  /**
  * map trimmed reads to assembled scaffolds
  */

  tag { sample_name }

  cpus 8
  memory '15 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  tuple path(trimmed_fq1), path(trimmed_fq2)
  path(scaffold_bt2index)

  output:
  path("${sample_name}.sorted.bam"), emit: bam
  path("${sample_name}_alnStats.txt"), emit: map_stats

  script:
  bam = "${sample_name}.sorted.bam"
  bai = "${sample_name}.bam.bai"
  stats = "${sample_name}_alnStats.txt"
  bt2_index = "${scaffold_bt2index}/${sample_name}"
  """
  echo MAPPING TO $bt2_index
  bowtie2 --very-sensitive -p ${task.cpus} -x ${scaffold_bt2index}/${sample_name} -1 $fq1 -2 $fq2 2> ${sample_name}_alnStats.txt | samtools view -h - | samtools sort - -o ${bam}
  """

  stub:
  bam = "${sample_name}.sorted.bam"
  stats_txt = "${sample_name}_alnStats.txt"

  """
  touch ${bam}
  touch ${stats_txt}
  """
}

process STRViper {
  /**
  * Process STRs discovered by Tandem Repeats Finder using STRViper
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/STR", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(bam)
  path(trf_dat)
  path(scaffolds)

  output:
  path("${sample_name}.vcf"), emit: str_vcf

  script:
  """
  jsat.str parseTRF --input ${trf_dat} --output ${sample_name}.trf.str --format str
  jsat.str sam2fragment --input ${bam} --output ${sample_name}.fragment
  jsat.str sortFragment --input ${sample_name}.fragment --output ${sample_name}.sorted.fragment
  jsat.str fragment2var --trfFile ${sample_name}.trf.str --output ${sample_name}.strv ${sample_name}.sorted.fragment
  jsat.str strv2vcf --input ${sample_name}.strv --output ${sample_name}.vcf --reference ${scaffolds}
  """

  stub:
  str_vcf = "${sample_name}.vcf"
  """
  touch ${str_vcf}
  """
}

process SNPEff {
  /**
  * Process STRs discovered by Tandem Repeats Finder using SNPEff
  */

  tag { sample_name }

  publishDir "${params.output_dir}/$sample_name/VariantAnalysis/STR", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  tuple val(sample_name), path(fq1), path(fq2)
  path(str_vcf)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  path("${sample_name}.ann.vcf"), emit: ann_vcf

  script:
  """
  java -jar \$SNPEFF_DIR/snpEff.jar -c \$SNPEFF_DIR/snpEff.config ${sample_name} ${str_vcf} > ${sample_name}.ann.vcf
  """

  stub:
  ann_vcf = "${sample_name}.ann.vcf"
  """
  touch ${ann_vcf}
  """
}
