process IQTREE {
    tag "${prefix}"
    label 'process_high'
    
    input:
    tuple val(meta), path(aln), path(const_sites), val(count)

    output:
    tuple val(meta), path("*.nwk"), emit: tree
    path 'versions.yml',            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args         = task.ext.args ?: ''
    prefix       = "${meta.species}-${meta.subtype}"
    uniq_seq = (count?.isInteger() ? count.toInteger() : 0)
    bootstrap    = uniq_seq > 4 ? '-B 1000' : ''
    """
    # run IQTREE3
    iqtree3 \\
        -s ${aln} \\
        -fconst \$(cat ${const_sites}) \\
        -T ${task.cpus} \\
        ${args} \\
        ${bootstrap}

    mv *.treefile ${prefix}.nwk

    #### VERSION INFO ####
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree2: \$(iqtree2 --version | head -n 1 | cut -f 4  -d ' ')
    END_VERSIONS
    """
}
