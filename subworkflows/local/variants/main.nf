//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ALIGN_READS        } from '../../../modules/local/align_reads/main'
include { SAMTOOLS_MPILEUP   } from '../../../modules/nf-core/samtools/mpileup/main'
include { LOWSITES           } from '../../../modules/local/lowsites/main'
include { FREEBAYES          } from '../../../modules/local/freebayes/main'
include { FILTER_VCF         } from '../../../modules/local/filter_vcf/main'
include { CREATE_MASK        } from '../../../modules/local/create_mask/main'
include { BCFTOOLS_CONSENSUS } from '../../../modules/local/bcftools/consensus/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL VARIANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow VARIANTS {

    take:
    // ch_samplesheet: [ val(meta), path(reads) ]
    ch_samplesheet
    // ch_refs: one record per reference, e.g. PREPARE.out.refs
    //          (channel.fromList(ReferenceManifest.load(...).records))
    ch_refs
    // make_consensus: boolean flag controlling consensus generation
    make_consensus

    main:

    // Collectors
    ch_versions = channel.empty()
    ch_aln      = channel.empty()

    // -------------------------------------------------------------------------
    // Samples with a reference supplied in the samplesheet keep it. All others
    // get the reference (and ploidy) whose species and subtype both match.
    // Samples with no matching reference are dropped, as before.
    // -------------------------------------------------------------------------
    ch_samplesheet
        .branch { meta, reads ->
            supplied: Utils.hasValue(meta.reference)
            lookup:   true
        }
        .set { ch_input }

    ch_input.lookup
        .combine( ch_refs.toList().map { refs -> [ refs ] } )
        .map { meta, reads, refs ->
            [ meta, reads, refs.find { ref -> matchesName(ref.species, meta.species) && matchesName(ref.subtype, meta.subtype) } ]
        }
        .filter { meta, reads, ref -> ref }
        .map { meta, reads, ref -> [ meta + [ reference: ref.reference, ploidy: ref.ploidy ], reads ] }
        .mix( ch_input.supplied )
        .set { ch_samplesheet }

    // -------------------------------------------------------------------------
    // ALIGN_READS
    // -------------------------------------------------------------------------
    ALIGN_READS(
        ch_samplesheet.map{ meta, reads -> [meta, reads, meta.reference] }
    )
    ch_versions = ch_versions.mix(ALIGN_READS.out.versions.first())
    ALIGN_READS.out.mapped.set { ch_mapped }

    // -------------------------------------------------------------------------
    // FREEBAYES (variant calling)
    // -------------------------------------------------------------------------
    FREEBAYES(
        ch_mapped.map{ meta, bam, bai -> [meta, bam, bai, meta.reference] }
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())

    // -------------------------------------------------------------------------
    // FILTER_VCF (post-calling filtering)
    // -------------------------------------------------------------------------
    FILTER_VCF(
        FREEBAYES.out.vcf
    )
    ch_versions = ch_versions.mix(FILTER_VCF.out.versions.first())

    // -------------------------------------------------------------------------
    // Optional: consensus generation + (optional) push to DB
    // When enabled, compute mpileup/low-sites for masking and build consensus.
    // -------------------------------------------------------------------------
    ch_depth = ch_samplesheet.map{[it[0], []]}
    if (make_consensus) {

        // SAMTOOLS mpileup (depth & pileup summaries)
        // mpileup inputs: [ meta, bam, opts ]; opts left as [] to preserve behavior
        SAMTOOLS_MPILEUP(
            ch_mapped.map { meta, bam, bai -> [ meta, bam, [] ] },
            [ null, [] ]
        )
        ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

        // LOWSITES (derive low-coverage sites and summary)
        LOWSITES(
            SAMTOOLS_MPILEUP.out.mpileup
        )
        ch_versions = ch_versions.mix(LOWSITES.out.versions.first())
        ch_depth    = LOWSITES.out.summary

        CREATE_MASK (
            LOWSITES.out.bed
                .join(FILTER_VCF.out.bed)
        )

        // Build consensus using filtered SNVs + reference + low-site mask
        BCFTOOLS_CONSENSUS(
            FILTER_VCF
                .out
                .snvs
                .map{ meta, vcf, tbi -> [meta, vcf, tbi, meta.reference] }
                .join(CREATE_MASK.out.bed)

        )
        ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())
        ch_aln      = BCFTOOLS_CONSENSUS.out.fasta

        // Optional: copy consensus into a DB layout by species/subtype
        if (params.push) {
            BCFTOOLS_CONSENSUS.out.fasta
                .subscribe { meta, fa ->
                    fa.copyTo(
                        file(params.db)
                            .resolve(meta.species)
                            .resolve(meta.subtype)
                            .resolve(fa.name)
                    )
                }
        }
    }

    emit:
    samplesheet = ch_samplesheet
    bam         = ch_mapped
    vcf         = FILTER_VCF.out.filt
    aln         = ch_aln
    depth       = ch_depth
    versions    = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Whether a sample value (a name or list of names) matches any of a reference
// record's names, compared after Utils.sanitize(). An unset value never matches.
//
def matchesName(names, value) {
    def wanted = (value instanceof Collection ? value : [ value ])
        .findAll { v -> Utils.hasValue(v) }
        .collect { v -> Utils.sanitize(v.toString()) }
    return wanted && names.any { n -> Utils.sanitize(n.toString()) in wanted }
}
