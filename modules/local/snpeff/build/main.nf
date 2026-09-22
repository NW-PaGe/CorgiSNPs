process SNPEFF_BUILD {
    tag "$meta.id"
    label 'process_medium'

    input:
    // fasta and gff may be gzipped or plain; they are staged under inputs/
    tuple val(meta), path(fasta, stageAs: 'inputs/*'), path(gff, stageAs: 'inputs/*')

    output:
    tuple val(meta), path("*.snpeff_db.tar.gz"), emit: db
    path "versions.yml"                        , emit: versions

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
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    # Database layout: snpeff_db/data/<genome>/{sequences.fa,genes.gff}
    # gzip -cdf decompresses gzipped files and passes plain files through
    mkdir -p snpeff_db/data/${meta.species}
    gzip -cdf ${fasta} > snpeff_db/data/${meta.species}/sequences.fa
    gzip -cdf ${gff}   > snpeff_db/data/${meta.species}/genes.gff

    cat <<-END_CONFIG > snpeff_db/snpEff.config
    data.dir = ./data/
    ${meta.species}.genome : ${meta.species}
    END_CONFIG

    snpEff \\
        -Xmx${avail_mem}M \\
        build \\
        -gff3 \\
        -noCheckCds \\
        -noCheckProtein \\
        -c snpeff_db/snpEff.config \\
        -dataDir \$PWD/snpeff_db/data \\
        ${args} \\
        ${meta.species}

    tar -czf ${prefix}.snpeff_db.tar.gz snpeff_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
