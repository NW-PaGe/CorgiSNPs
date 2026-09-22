process EXTRACT_REGIONS {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(bam), path(ref), path(genome), path(gff), val(genes)

    output:
    tuple val(meta), path("*.fastq.gz"), path("genome_masked.fasta.gz"), emit: results
    path "versions.yml",                                                 emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def gene_args = genes.collect { g -> "'${g}'" }.join(' ')
    tool = 'coordcutter'
    """
    ${tool} \\
        --bam ${bam} \\
        --ref ${ref} \\
        --genome ${genome} \\
        --gff ${gff} \\
        --genes ${gene_args}

    gzip genome_masked.fasta
    
    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${tool}: "\$(${tool} --version 2>&1 | tr -d '\\r')"
    END_VERSIONS
    """
}