process GLNEXUS {
    tag "${meta.id}"
    label 'process_high'

    // GLnexus does not support Conda: https://github.com/dnanexus-rnd/GLnexus/issues/74
    container "quay.io/mlin/glnexus:v1.3.1"

    input:
    tuple val(meta), path(gvcfs), path(tbis), path(intervals)

    output:
    tuple val(meta), path("${prefix}.bcf"), emit: bcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GLNEXUS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def bed    = intervals ? "--bed ${intervals}" : ''
    // GLnexus's internal cache is capped below the task memory so the
    // process itself (and htslib decompression buffers) have headroom
    def mem_gbytes = Math.max(1, (task.memory.toGiga() * 0.9) as int)
    """
    ulimit -n 65536

    # Fresh scratch DB dir per attempt: glnexus_cli refuses to start if
    # GLnexus.DB already exists (relevant on task.attempt retries)
    rm -rf GLnexus.DB

    # Write the gVCF list to a file for reproducibility/debugging. At ~1000
    # samples the resulting argument list (~50KB of paths) stays well under
    # the typical 2MB ARG_MAX, so it is passed directly via \$(cat ...).
    printf '%s\\n' ${gvcfs} > gvcf.list

    glnexus_cli \\
        --config DeepVariant \\
        --threads ${task.cpus} \\
        --mem-gbytes ${mem_gbytes} \\
        ${bed} \\
        ${args} \\
        \$(cat gvcf.list) > ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --version 2>&1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: \$(glnexus_cli --version 2>&1 | sed 's/^.*release v//; s/ .*\$//')
    END_VERSIONS
    """
}
