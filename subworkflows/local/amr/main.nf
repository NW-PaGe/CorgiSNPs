//
// Subworkflow for AMR analysis in the mycosnp pipeline
//

include { SNPEFF_BUILD    } from '../../../modules/local/snpeff/build/main'
include { EXTRACT_REGIONS } from '../../../modules/local/extract_regions/main'
include { VARIANTS        } from '../../../subworkflows/local/variants'
include { SNPEFF          } from '../../../modules/local/snpeff/main'
include { SNPEFF_PARSE    } from '../../../modules/local/snpeff_parse/main'


/*
    Each sample is routed by the manifest record for its reference:

    direct   The sample's reference has 'amr' (and so 'annotation'). Its VCF is
             annotated with a snpEff database built from that reference.
    extract  The sample's reference has no 'amr', or isn't in the manifest.
             Target regions are extracted against the species' primary
             reference, variants are re-called, and the new VCF is annotated
             with the primary's snpEff database.
    none     The species has no reference with 'amr'. AMR is skipped for the
             sample, with a warning.

    A snpEff database is built once per reference that is used, and the
    reference's 'amr' block is written out as the JSON for SNPEFF_PARSE.
*/
workflow AMR {

    take:
    // ch_samplesheet: [ meta, reads ]
    ch_samplesheet
    // ch_bam : [ meta, bam, bai ]
    ch_bam
    // ch_vcf : [ meta, vcf, csi ]
    ch_vcf
    // ch_refs: one record per reference, e.g. PREPARE.out.refs
    //          (channel.fromList(ReferenceManifest.load(...).records))
    ch_refs

    main:
    // Collector for version files
    ch_versions = channel.empty()

    // ---------------------------------------------------------------------
    // Route each sample: [ meta, amr_reference_record, is_direct ]
    // ---------------------------------------------------------------------
    ch_samplesheet
        .combine( ch_refs.toList().map { refs -> [ refs ] } )
        .map { meta, reads, refs ->
            def own    = refs.find { ref -> sameReference(ref, meta) }
            def direct = own != null && ReferenceManifest.hasAmr(own)
            def target = direct ? own : refs.find { ref -> ref.primary && matchesSpecies(ref, meta.species) }
            [ meta, target, direct ]
        }
        .branch { meta, ref, direct ->
            direct:  direct
            extract: ref
            none:    true
        }
        .set { ch_route }

    ch_route.none
        .map { meta, ref, direct -> "${meta.id} (${meta.species})".toString() }
        .collect()
        .subscribe { samples ->
            log.warn "AMR skipped for ${samples.size()} sample(s) with no reference AMR targets for their species: ${samples.join(', ')}"
        }

    // References used for AMR, once each
    ch_route.direct
        .mix( ch_route.extract )
        .map { meta, ref, direct -> ref }
        .unique { ref -> refKey(ref) }
        .set { ch_amr_refs }

    // ---------------------------------------------------------------------
    // Build a snpEff database per reference: [ [id: ref_key], db.tar.gz ]
    // ---------------------------------------------------------------------
    SNPEFF_BUILD(
        ch_amr_refs.map { ref ->
            [ [ id: refKey(ref), species: Utils.sanitize(ref.species[0]) ], ref.reference, ref.annotation ]
        }
    )
    ch_versions = ch_versions.mix(SNPEFF_BUILD.out.versions.first())

    // AMR targets per reference, as JSON for SNPEFF_PARSE: [ ref_key, json ]
    ch_amr_refs
        .collectFile { ref ->
            [ "${refKey(ref)}.amr.json", groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(ref.amr)) + '\n' ]
        }
        .map { json -> [ json.name - '.amr.json', json ] }
        .set { ch_amr_json }

    // ---------------------------------------------------------------------
    // Direct: annotate the original VCF
    // ---------------------------------------------------------------------
    ch_route.direct
        .map  { meta, ref, direct -> [ meta, refKey(ref) ] }
        .join ( ch_vcf.map { meta, vcf, csi -> [ meta, vcf ] } )
        .map  { meta, key, vcf -> [ meta + [ amr_ref: key ], vcf ] }
        .set  { ch_direct_vcf }

    // ---------------------------------------------------------------------
    // Extract: pull target regions against the primary reference and
    // re-call variants
    // ---------------------------------------------------------------------
    EXTRACT_REGIONS(
        ch_route.extract
            .map  { meta, ref, direct -> [ meta, ref ] }
            .join ( ch_bam.map { meta, bam, bai -> [ meta, bam ] } )
            .map  { meta, ref, bam ->
                [ meta + [ amr_ref: refKey(ref) ], bam, meta.reference, ref.reference, ref.annotation ]
            }
    )
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions.first())

    VARIANTS(
        EXTRACT_REGIONS.out.results.map { meta, reads, ref -> [ meta + [ reference: ref ], reads ] },
        ch_refs,
        false
    )
    ch_versions = ch_versions.mix(VARIANTS.out.versions)

    // ---------------------------------------------------------------------
    // SnpEff: annotate each VCF with its reference's database
    // ---------------------------------------------------------------------
    SNPEFF(
        ch_direct_vcf
            .mix( VARIANTS.out.vcf.map { meta, vcf, csi -> [ meta, vcf ] } )
            .map { meta, vcf -> [ meta.amr_ref, meta, vcf ] }
            .combine( SNPEFF_BUILD.out.db.map { meta, db -> [ meta.id, db ] }, by: 0 )
            .map { key, meta, vcf, db -> [ meta, vcf, db ] }
    )
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())

    // Parse SnpEff results with the reference's AMR targets
    SNPEFF_PARSE(
        SNPEFF.out.vcf
            .map { meta, vcf -> [ meta.amr_ref, meta, vcf ] }
            .combine( ch_amr_json, by: 0 )
            .map { key, meta, vcf, json -> [ meta, vcf, json ] }
    )
    ch_versions = ch_versions.mix(SNPEFF_PARSE.out.versions.first())

    emit:
    summary   = SNPEFF_PARSE.out.target
    snpeff_db = SNPEFF_BUILD.out.db   // [ [id: ref_key, species], snpeff_db.tar.gz ]
    versions  = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Unique key for a reference record, e.g. 'candidozyma_auris_clade_i'
//
def refKey(ref) {
    return "${ref.name}_${Utils.sanitize(ref.subtype[0].toString())}".toString()
}

//
// Whether a record is the reference the sample was aligned to
//
def sameReference(ref, meta) {
    return ref.reference && Utils.hasValue(meta.reference) && ref.reference.toString() == meta.reference.toString()
}

//
// Whether a sample's species matches any of a record's species names
//
def matchesSpecies(ref, species) {
    def wanted = (species instanceof Collection ? species : [ species ])
        .findAll { v -> Utils.hasValue(v) }
        .collect { v -> Utils.sanitize(v.toString()) }
    return wanted && ref.species.any { n -> Utils.sanitize(n.toString()) in wanted }
}
