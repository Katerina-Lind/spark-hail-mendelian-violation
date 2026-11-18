nextflow.enable.dsl = 2

params.chunks = null
params.pedigree = null
params.lcr_bed = null
params.results = null
params.hail_tmp = null
params.hail_container = '249640142434.dkr.ecr.us-east-1.amazonaws.com/batch/hail-mendelian:0.2.138-s3-v1'

process NATIVE_DE_NOVO {
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
    """
    export TMPDIR=/tmp
    export SPARK_LOCAL_DIRS=/tmp
    export PYSPARK_SUBMIT_ARGS='--driver-memory 40g --executor-memory 40g pyspark-shell'
    python '${runner}' \\
      --chunk-id '${chunk_id}' \\
      --matrix '${matrix}' \\
      --prep-manifest '${prep_manifest}' \\
      --pedigree '${params.pedigree}' \\
      --lcr-bed '${params.lcr_bed}' \\
      --tmp-base '${params.hail_tmp}' \\
      --out-prefix '${params.results}/${chunk_id}'
    """
}

workflow {
    ['chunks', 'pedigree', 'lcr_bed', 'results', 'hail_tmp'].each { p ->
        if (!params[p]) error "Missing --${p}"
    }
    chunks_ch = Channel.fromPath(params.chunks, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.chunk_id, row.matrix, row.manifest) }
    script_ch = Channel.value(file("${projectDir}/scripts/03_hail_de_novo.py", checkIfExists: true))
    NATIVE_DE_NOVO(chunks_ch, script_ch)
}
