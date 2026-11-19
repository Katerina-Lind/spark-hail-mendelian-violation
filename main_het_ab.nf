nextflow.enable.dsl = 2

params.chunks = null
params.samples = 'SDSM-V8W,SDSM-2FNA,SDSM-213R,SDSM-213N,SDSM-2116'
params.results = null
params.hail_tmp = null
params.pass_only = false
params.hail_container = '249640142434.dkr.ecr.us-east-1.amazonaws.com/batch/hail-mendelian:0.2.138-s3-v1'

process EXPORT_HET_AB {
    tag "$chunk_id"
    cpus 8
    memory '64 GB'
    container "${params.hail_container}"
    errorStrategy 'terminate'
    input:
    tuple val(chunk_id), val(matrix), val(prep_manifest)
    path runner
    output:
    val(chunk_id)
    script:
    def passArg = params.pass_only ? '--pass-only' : ''
    """
    export TMPDIR=/tmp
    export SPARK_LOCAL_DIRS=/tmp
    export PYSPARK_SUBMIT_ARGS='--driver-memory 40g --executor-memory 40g pyspark-shell'
    python '${runner}' \\
      --chunk-id '${chunk_id}' \\
      --matrix '${matrix}' \\
      --prep-manifest '${prep_manifest}' \\
      --samples '${params.samples}' \\
      --out-prefix '${params.results}/${chunk_id}' \\
      --tmp-base '${params.hail_tmp}' ${passArg}
    """
}

workflow {
    ['chunks','samples','results','hail_tmp'].each { p ->
        if (!params[p]) error "Missing --${p}"
    }
    chunks_ch = Channel.fromPath(params.chunks, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.chunk_id,row.matrix,row.manifest) }
    runner_ch = Channel.value(file("${projectDir}/scripts/04_export_het_snp_ab.py", checkIfExists: true))
    EXPORT_HET_AB(chunks_ch, runner_ch)
}
