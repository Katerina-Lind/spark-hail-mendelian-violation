nextflow.enable.dsl = 2

params.chunks = null
params.samples = 'SDSM-V8W,SDSM-V8Y,SDSM-V8V,SDSM-2FNA,SDSM-2GIZ,SDSM-2FNB'
params.results = 's3://ucsf-spark/01-sfari/analysis/trio-wgs/hail-mendelian/results/kinship_chr1_chr2_chr22_v1'
params.hail_tmp = 's3://ucsf-spark/01-sfari/analysis/trio-wgs/hail-mendelian/tmp'
params.hail_container = '249640142434.dkr.ecr.us-east-1.amazonaws.com/batch/hail-mendelian:0.2.138-s3-v1'

process EXTRACT_KINSHIP_MARKERS {
    tag "$chunk_id"
    cpus 8
    memory '64 GB'
    container "${params.hail_container}"
    errorStrategy 'terminate'
    input:
    tuple val(chromosome), val(chunk_id), val(matrix), val(prep_manifest)
    path runner
    output:
    val(chunk_id)
    script:
    """
    export TMPDIR=/tmp
    export SPARK_LOCAL_DIRS=/tmp
    export PYSPARK_SUBMIT_ARGS='--driver-memory 40g --executor-memory 40g pyspark-shell'
    python '${runner}' \\
      --chromosome '${chromosome}' \\
      --chunk-id '${chunk_id}' \\
      --matrix '${matrix}' \\
      --prep-manifest '${prep_manifest}' \\
      --samples '${params.samples}' \\
      --tmp-base '${params.hail_tmp}' \\
      --out-prefix '${params.results}/${chunk_id}'
    """
}

process SUMMARISE_KINSHIP {
    tag 'chr1_chr2_chr22'
    cpus 2
    memory '8 GB'
    container "${params.hail_container}"
    errorStrategy 'terminate'
    input:
    val completed_ids
    path runner
    output:
    val('kinship_complete')
    script:
    """
    python '${runner}' \\
      --manifest '${params.chunks}' \\
      --samples '${params.samples}' \\
      --results '${params.results}'
    """
}

workflow {
    ['chunks','samples','results','hail_tmp'].each { p ->
        if (!params[p]) error "Missing --${p}"
    }
    chunks_ch = Channel.fromPath(params.chunks, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.chromosome, row.chunk_id, row.matrix, row.manifest) }
    extract_script = Channel.value(file("${projectDir}/scripts/05_extract_kinship_markers.py", checkIfExists: true))
    summary_script = Channel.value(file("${projectDir}/scripts/06_summarise_kinship.py", checkIfExists: true))
    extracted = EXTRACT_KINSHIP_MARKERS(chunks_ch, extract_script)
    SUMMARISE_KINSHIP(extracted.collect(), summary_script)
}
