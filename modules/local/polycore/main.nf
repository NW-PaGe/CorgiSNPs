process POLYCORE {
    tag "${prefix}"
    label 'process_high'
        
    input:
    tuple val(meta), path(ref), path(assemblies)

    output:
    tuple val(meta), path("core.aln.gz"),       emit: snps
    tuple val(meta), path("core.full.aln.gz"),  emit: full
    tuple val(meta), path("summary.csv"),       emit: csv
    tuple val(meta), path("full.csv.gz"),       emit: full_csv
    tuple val(meta), path("fconst.txt"),        emit: fconst
    tuple val(meta), path("dist_long.csv.gz"), emit: dist_long
    tuple val(meta), path("dist_wide.csv.gz"),  emit: dist_wide
    tuple val(meta), path("*.html"),            emit: plot, optional: true
    tuple val(meta), env("UNIQ_SEQ"),           emit: uniq_seq
    path "versions.yml",                        emit: versions


    script:
    def args = task.ext.args ?: ''
    prefix = "${meta.species}-${meta.subtype}"
    def origName = ref.getName()
    def stem = origName.replaceAll(/\.gz$/, '').replaceAll(/\.(fna|fa|fasta|fas)$/, '')
    def suffix = origName.substring(stem.length())
    def prefixedName = stem.startsWith('Reference_') ? origName : "Reference_${stem}${suffix}"
    """
    if [ "${origName}" != "${prefixedName}" ]; then
        ln -s ${ref} ${prefixedName}
    fi

    polycore \\
        ${prefixedName} ${assemblies} \\
        --min-gf ${params.min_genome_fraction} \\
        --min-cf ${params.min_core_fraction} \\
        --ploidy ${meta.ploidy} \\
        ${args}

    UNIQ_SEQ=\$(cat core.aln | grep -v '>' | sort | uniq | wc -l)

    gzip *.aln dist_*.csv full.csv
    
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polycore: "\$(polycore --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}
