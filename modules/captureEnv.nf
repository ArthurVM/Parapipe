process captureEnv {
  /**
  * Capture the run environment
  */

  publishDir "${params.output_dir}", mode: 'copy', overwrite: 'true'

  memory '5 GB'

  input:
    val(input_dir)
    val(output_dir)
    val(ref)
    val(yaml)
    val(db)
    val(read_n_threshold)
    val(mincov)
    val(missing)
    val(maf)
    val(mac)


  output:
    path("env.json"), emit: env_json

  script:
    scripts = "${baseDir}/bin"
    """
    python3 ${scripts}/get_environment.py ${input_dir} ${output_dir} ${ref} ${yaml} ${db} ${baseDir}/singularity/ ${read_n_threshold} ${mincov} ${missing} ${maf} ${mac}
    """
    
  stub:
    """
    touch env.json
    """
}
