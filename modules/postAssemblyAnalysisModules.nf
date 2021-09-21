// modules for post assembly genome analysis

process makeChromosomeFastas {
  /**
  * generates fastas of individual chromosomes from a set of assembly fastas scaffolded using ABACAS
  */

  memory '5 GB'

  input:
  val(assemblies)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  // val("\$PWD"), emit: wd
  path("*fasta"), emit: chr_multifastas

  script:
  """
  echo ${assemblies}
  python3 ${baseDir}/scripts/makeChromosomeFastas.py ${assemblies}
  """

  stub:
  """
  touch {chr0,chr1,chr2,chr3,chr4}.fasta
  """
}

process chromosomeAlignment {
  /**
  * Align chromosomes using preogressiveMauve
  */

  publishDir "${params.output_dir}/ChrAln", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
  path(chr_multifastas)

  output:
  path("*xmfa"), emit: xmfa

  script:
  """
  chr_mfs="${chr_multifastas}"
  IFS=', ' read -r -a chr_mf_array <<< "\${chr_mfs}"
  for chr in \${chr_mf_array[@]}
  do
    echo \$chr
    fname="\${chr##*/}"
    chrID="\${fname%.*}"
    progressiveMauve --output=\${chrID}.xmfa \$chr
  done
  """

  stub:
  """
  chr_mfs="${chr_multifastas}"
  IFS=', ' read -r -a chr_mf_array <<< "\${chr_mfs}"
  for chr in \${chr_mf_array[@]}
  do
    fname="\${chr##*/}"
    chrID="\${fname%.*}"
    echo "progressiveMauve --output=\${chrID}.xmfa \$chr"
    touch \${chrID}.xmfa
  done
  """
}
