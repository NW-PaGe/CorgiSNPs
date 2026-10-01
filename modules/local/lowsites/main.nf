process LOWSITES {
    tag "$meta.id"
    label 'process_low'
        
    input:
    tuple val(meta), path(pileup)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.csv"), emit: summary
    path "versions.yml",            emit: versions

    script:
    def args = task.ext.args ?: '' 
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cfg    = meta.settings ?: params
    tool = 'lowsites'
    """
    ${tool} \\
        -t ${cfg.min_base_quality} \\
        -d ${cfg.min_base_depth} \\
        -o ${prefix}.mask.bed \\
        ${pileup}

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
