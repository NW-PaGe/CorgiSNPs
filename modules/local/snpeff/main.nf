process SNPEFF {
    tag "$meta.id"
    label 'process_medium'

    input:
    // db: snpeff_db.tar.gz from SNPEFF_BUILD
    tuple val(meta), path(vcf), path(db)

    output:
    tuple val(meta), path("*.ann.vcf"), emit: vcf
    tuple val(meta), path("*.csv")    , emit: report
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tar -xzf ${db}

    # genome ID is the '<id>.genome : <name>' entry in the database's config
    genome=\$(sed -nE 's/^([^ ]+)\\.genome *:.*/\\1/p' snpeff_db/snpEff.config | head -n 1)

    snpEff \\
        -Xmx${avail_mem}M \\
        \$genome \\
        -c snpeff_db/snpEff.config \\
        -dataDir \$PWD/snpeff_db/data \\
        -csvStats ${prefix}.csv \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.ann.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ann.vcf
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
