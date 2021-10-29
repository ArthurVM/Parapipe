// modules for post assembly genome analysis

process makeChromosomeFastas {
  /**
  * generates fastas of individual chromosomes from a set of assembly fastas scaffolded using ABACAS
  */

  memory '5 GB'

  input:
  path(assemblies)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:
  // val("\$PWD"), emit: wd
  path("*fasta"), emit: chr_multifastas

  script:
  """
  echo ${assemblies}
  python3 ${baseDir}/scripts/makeChromosomeFastas.py ${assemblies} ${fasta}
  """

  stub:
  """
  python3 ${baseDir}/scripts/makeChromosomeFastas.py ${assemblies} ${fasta}
  """
}

process chromosomeAlignment {
  /**
  * Align chromosomes using preogressiveMauve
  */

  publishDir "${params.output_dir}/ChrAln", mode: 'copy', pattern: '*xmfa', overwrite: 'true'

  memory '15 GB'

  input:
  path(chr_multifastas)

  output:
  path("*xmfa"), emit: xmfa
  path("*.maln.map"), emit: maln_maps

  script:
  """
  chr_mfs="${chr_multifastas}"
  IFS=', ' read -r -a chr_mf_array <<< "\${chr_mfs}"
  for chr in \${chr_mf_array[@]}
  do
    echo \$chr
    fname="\${chr##*/}"
    chrID="\${fname%.*}"
    echo "progressiveMauve --output=\${chrID}.xmfa \$chr"
    progressiveMauve --output=\${chrID}.xmfa \$chr
    grep \">\" \$chr | cat -n | sed -e \"s/>//g\" | awk \'{{split(\$0, a ,\" \");print (a[1]\" \"a[2])}}\' > \${chrID}.maln.map
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
    grep \">\" \$chr | cat -n | sed -e \"s/>//g\" | awk \'{{split(\$0, a ,\" \");print (a[1]\" \"a[2])}}\' > \${chrID}.maln.map
  done
  """
}

process dNdS {
  /**
  * Aligns orthologs and calculates selection pressure for each gene
  */

  memory '10 GB'

  input:
  path(assemblies)
  path(annotations)
  tuple path(fasta), path(gff), path(cds), path(gaf)

  output:

  script:
  """
  python3 ${baseDir}/scripts/gff2protein.py -g ${annotations} -f ${assemblies} -id gene
  python3 ${baseDir}/scripts/alignOrthos.py ./*-cdna.fa
  """

  stub:
  """
  echo "python3 ${baseDir}/scripts/gff2protein.py -g ${annotations} -f ${assemblies} -id gene"
  echo "python3 ${baseDir}/scripts/alignOrthos.py ./*-cdna.fa"
  """
}
