process SUBTYPE {
    tag "$meta.id"
    label 'process_low'

    input:
    // signatures: the species' reference signatures from SUBTYPE_SKETCH, named by subtype
    // threshold:  ANI threshold for the species (fraction or percent); '' uses the script default
    tuple val(meta), path(contigs), path(signatures), val(threshold)

    output:
    tuple val(meta), path("*_subtype.csv"), emit: subtype
    path "versions.yml"                   , emit: versions

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def tool          = 'subtyper.py'
    def threshold_arg = threshold != null && threshold.toString().trim() ? "--threshold ${threshold}" : ''
    """
    ${tool} assign \\
        ${signatures} \\
        ${contigs} \\
        --sample ${meta.id} \\
        --out ${prefix}_subtype.csv \\
        ${threshold_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample,subtype\\n${meta.id},undefined\\n' > ${prefix}_subtype.csv
    printf '"${task.process}":\\n    subtyper.py: stub\\n' > versions.yml
    """
}
