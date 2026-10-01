process FILTER_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5acacb55c52bec97c61fd34ffa8721fce82ce823005793592e2a80bf71632cd0/data':
        'community.wave.seqera.io/library/bcftools:1.21--4335bec1d7b44d11' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tagged.vcf.gz"), path("*.tagged.vcf.gz.tbi"), emit: tagged
    tuple val(meta), path("*.filt.vcf.gz"), path("*.filt.vcf.gz.tbi"),     emit: filt
    tuple val(meta), path("*.snvs.vcf.gz"), path("*.snvs.vcf.gz.tbi"),     emit: snvs
    tuple val(meta), path("*.bed"),                                        emit: bed
    path  "versions.yml"                                             ,     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Thresholds from the sample's reference settings (subtype > species > run)
    def cfg    = meta.settings ?: params

    """
    bcftools index -t ${vcf}

    # 1) Compute per-sample VAF using ALL alleles (vector before split)
    bcftools +fill-tags ${vcf} -Ov -- -t FORMAT/VAF | \\
    # 2) Split multiallelics to one ALT per line
    bcftools norm -m -any -Ou | \\
    # 3) Fill allele count (AC) from genotype (GT)
    bcftools +fill-tags -Ou -- -t AC | \\
    # 4) Tag non-GT alleles
    bcftools filter -Ou -s non_gt_alt -e 'INFO/AC==0' | \\
    # 5) *** AF filter against ALL detected alleles ***
    bcftools filter -Ou -s low_af -e 'FMT/VAF < ${cfg.min_allele_fraction}' | \\
    # 6) Tag low allele depth
    bcftools filter -Ou -s low_depth -e 'FMT/AO < ${cfg.min_base_depth}' | \\
    # 7) Tag low allele base quality
    bcftools filter -Ou -s low_base_qual -e '(FMT/AO > 0 && (FMT/QA)/(FMT/AO) <= ${cfg.min_base_quality})' | \\
    # 8) Tag low allele mapping quality
    bcftools filter -Ou -s low_map_qual -e 'INFO/MQM < ${cfg.min_mapping_quality}' | \\
    # 9) Tag strand bias
    bcftools filter -Ou -s strand_bias -e '(INFO/SAP >= ${cfg.max_strand_bias}) || (INFO/SAF+INFO/SAR > 0 && ((INFO/SAF)/(INFO/SAF+INFO/SAR) < ${cfg.min_fwd_strand_fraction} ))' | \\
    # 10) Read position bias
    bcftools filter -Ou -s read_pos_bias -e '(INFO/RPR=0 || INFO/RPL=0) || (INFO/RPP > ${cfg.max_read_pos_bias})' | \\
    # 11) Save to VCF
    bcftools view -Oz -o ${prefix}.tagged.vcf.gz

    bcftools index -t ${prefix}.tagged.vcf.gz

    #-- FOR AMR VARIANT CALLING --#
    # Remove filtered sites
    bcftools view \\
        -i 'FILTER=="PASS"' \\
        -Oz \\
        -o ${prefix}.filt.vcf.gz \\
        ${prefix}.tagged.vcf.gz
    
    bcftools index -t ${prefix}.filt.vcf.gz

    #-- FOR CONSENSUS GENERATION --#
    # Filter to single variants only
    bcftools view \\
        -i 'TYPE=="snp"' \\
        -Oz \\
        -o ${prefix}.snvs.vcf.gz \\
        ${prefix}.filt.vcf.gz
    
    bcftools index -t ${prefix}.snvs.vcf.gz

    # Create mask for all sites that are not SNPs or did not pass a filter
    bcftools view \\
        -e 'FILTER="PASS" && TYPE="snp"' \\
        ${prefix}.tagged.vcf.gz | \\
        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' \\
        > ${prefix}.filt_vcf.bed

    # Create a BED of all PASS SNP positions (1 bp each after normalization)
    bcftools view \\
        -i 'FILTER="PASS" && TYPE="snp"' \\
        ${prefix}.tagged.vcf.gz | \\
        bcftools query -f '%CHROM\\t%POS0\\t%END\\n' \\
        > ${prefix}.pass_snp.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
