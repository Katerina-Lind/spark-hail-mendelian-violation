nextflow.enable.dsl = 2

params.chunks = null           // TSV: chunk_id,vcf,input_version,matrix,manifest
params.pedigree = null
params.lcr_bed = null
params.skip_lcr = false
params.results = null          // Hail-accessible S3 prefix, e.g. s3a://bucket/project/results
params.hail_tmp = null         // Hail-accessible shared scratch base
params.hail_container = null   // Pin digest of validated Hail container
params.emit_vcf = false

process PREPARE_MATRIX {
    tag { chunk_id }
    container params.hail_container
    cpus 8
    memory '64 GB'
    errorStrategy 'terminate'

    input:
    tuple val(chunk_id), val(vcf), val(input_version), val(matrix), val(manifest)
    path prep_script

    output:
    tuple val(chunk_id), val(matrix), val(manifest)

    script:
    """
export PYSPARK_SUBMIT_ARGS='--driver-memory 40g --executor-memory 40g pyspark-shell'
    python ${prep_script} \\
      --vcf '${vcf}' \\
      --matrix '${matrix}' \\
      --manifest '${manifest}' \\
      --input-version '${input_version}' \\
      --local-cores ${task.cpus} \\
      --tmp-dir '${params.hail_tmp}'
    """
}

process MENDELIAN_SCREEN {
    tag { chunk_id }
    container params.hail_container
    cpus 8
    memory '64 GB'
    errorStrategy 'terminate'

    input:
    tuple val(chunk_id), val(matrix), val(manifest)
    path pedigree
    path lcr_bed
    path analysis_script

    output:
    path 'completed.txt'

    script:
    def lcr_arg = params.skip_lcr ? '--skip-lcr' : "--lcr-bed '${lcr_bed}'"
    def vcf_arg = params.emit_vcf ? '' : '--no-vcf'
    """
    export PYSPARK_SUBMIT_ARGS='--driver-memory 40g --executor-memory 40g pyspark-shell'
    python ${analysis_script} \\
      --matrix '${matrix}' \\
      --prep-manifest '${manifest}' \\
      --pedigree '${pedigree}' \\
      --out-prefix '${params.results}/${chunk_id}' \\
      --local-cores ${task.cpus} \\
      --tmp-dir '${params.hail_tmp}' \\
      ${lcr_arg} ${vcf_arg}
    echo "${chunk_id}\t${params.results}/${chunk_id}.completed.json" > completed.txt
    """
}

workflow {
    for (p in ['chunks', 'pedigree', 'results', 'hail_tmp', 'hail_container']) {
        if (!params[p]) error "Missing --${p.replace('_','-')}"
    }
    if (!params.skip_lcr && !params.lcr_bed) error 'Supply --lcr_bed or --skip_lcr true'

    def csv = Channel.fromPath(params.chunks, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { r ->
            tuple(r.chunk_id, r.vcf, r.input_version, r.matrix, r.manifest)
        }

    prepared = PREPARE_MATRIX(csv, file("${projectDir}/scripts/01_prepare_matrix.py"))
    MENDELIAN_SCREEN(
        prepared,
        file(params.pedigree, checkIfExists: true),
        params.skip_lcr ? file("${projectDir}/README.md") : file(params.lcr_bed, checkIfExists: true),
        file("${projectDir}/scripts/02_mendelian_screen.py")
    )
}
