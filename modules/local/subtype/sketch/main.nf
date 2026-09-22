process SUBTYPE_SKETCH {
    tag "$meta.id"
    label 'process_low'

    input:
    // One item per species: the species' reference records and assemblies.
    // ref_meta and references are matched by file name (ref_meta.assembly).
    tuple val(meta), val(ref_meta), path(references, stageAs: 'staged/*')

    output:
    tuple val(meta), path("*.sig.zip"), emit: signatures
    path "versions.yml"               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tool   = 'subtyper.py'

    // subtypes.tsv maps each staged assembly to its subtype name. Built as a
    // single line (\t and \n expanded by printf %b) so the script block keeps
    // a consistent indent; multi-line interpolation breaks Nextflow's indent
    // stripping and with it the END_VERSIONS heredoc below.
    def table = ( [ 'subtype\\tpath' ] + ref_meta.collect { ref ->
        def name = ref.subtype instanceof List ? ref.subtype[0] : ref.subtype
        "${name}\\tstaged/${ref.assembly}"
    } ).join('\\n')

    """
    printf '%b\\n' '${table}' > subtypes.tsv

    ${tool} sketch \\
        subtypes.tsv \\
        ${references} \\
        --prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sig.zip
    printf '"${task.process}":\\n    subtyper.py: stub\\n' > versions.yml
    """
}
