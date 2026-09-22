//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SHOVILL        } from '../../../modules/local/shovill/main'
include { GAMBIT_QUERY   } from '../../../modules/local/gambit/main'
include { SUBTYPE_SKETCH } from '../../../modules/local/subtype/sketch/main'
include { SUBTYPE        } from '../../../modules/local/subtype/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CLASSIFY SAMPLES AND MATCH THEM TO REFERENCES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    1. Every sample is assembled with SHOVILL.
    2. Samples without meta.species go through GAMBIT to get one.
    3. Every sample is joined to the references for its species.
    4. For each species with samples to subtype, its references are sketched
       once into signatures named by subtype (SUBTYPE_SKETCH).
    5. Samples without meta.subtype are subtyped against their species'
       signatures (SUBTYPE).
    6. Each sample is matched to the reference for its species and subtype.

    Any sample that cannot be matched to a reference stops the pipeline, with
    every such sample and the reason listed in one error.
*/
workflow CLASSIFY {

    take:
    // ch_samplesheet: [ val(meta), path(reads) ]; meta.species / meta.subtype may already be set
    ch_samplesheet
    // ch_refs: one record per reference, e.g. channel.fromList(ReferenceManifest.load(...).records)
    ch_refs

    main:

    // Collector for version files
    ch_versions = channel.empty()

    // -------------------------------------------------------------------------
    // MODULE: de novo assembly (all samples)
    // -------------------------------------------------------------------------
    SHOVILL(
        ch_samplesheet
    )

    ch_samplesheet
        .join(SHOVILL.out.contigs)
        .branch { meta, reads, contigs ->
            has_species:   Utils.hasValue(meta.species)
            needs_species: true
        }
        .set { ch_assembled }

    // -------------------------------------------------------------------------
    // MODULE: Species assignment via GAMBIT (only samples without a species)
    // -------------------------------------------------------------------------
    GAMBIT_QUERY(
        ch_assembled.needs_species.map { meta, reads, contigs -> [ meta, contigs ] },
        params.gambit_db,
        params.gambit_h5_dir
    )
    ch_versions = ch_versions.mix(GAMBIT_QUERY.out.versions.first())

    // -------------------------------------------------------------------------
    // Add the GAMBIT species, combine with the samples that already had one,
    // and join each sample to the references for its species:
    // [ meta, reads, contigs, [ matching reference records ] ]
    // -------------------------------------------------------------------------
    ch_assembled.needs_species
        .join( GAMBIT_QUERY.out.taxa.map { meta, csv -> [ meta, gambitSpecies(csv) ] } )
        .map { meta, reads, contigs, species -> [ meta + [ species: species ], reads, contigs ] }
        .mix( ch_assembled.has_species )
        .combine( ch_refs.toList().map { refs -> [ refs ] } )
        .map { meta, reads, contigs, refs ->
            [ meta, reads, contigs, refs.findAll { ref -> matchesSpecies(ref, meta.species) } ]
        }
        .branch { meta, reads, contigs, refs ->
            no_species:    !Utils.hasValue(meta.species)
            no_reference:  !refs
            has_subtype:   Utils.hasValue(meta.subtype)
            needs_subtype: true
        }
        .set { ch_species }

    // -------------------------------------------------------------------------
    // MODULE: Sketch the references once per species that needs subtyping.
    // Input: [ [id: species_key], [ ref_meta ], [ reference assemblies ] ]
    // -------------------------------------------------------------------------
    SUBTYPE_SKETCH(
        ch_species.needs_subtype
            .map { meta, reads, contigs, refs -> [ speciesKey(refs), refs ] }
            .unique { it[0] }
            .map { key, refs ->
                [ [ id: key ], refs.collect { ref -> referenceInfo(ref) }, refs.collect { ref -> ref.reference } ]
            }
    )
    ch_versions = ch_versions.mix(SUBTYPE_SKETCH.out.versions.first())

    // -------------------------------------------------------------------------
    // MODULE: Subtyping (only samples without a subtype). Each sample is paired
    // with its species' signatures by species key, and gets the species' ANI
    // threshold from the manifest.
    // Input: [ meta, contigs, signatures, threshold ]
    // -------------------------------------------------------------------------
    SUBTYPE(
        ch_species.needs_subtype
            .map { meta, reads, contigs, refs -> [ speciesKey(refs), meta, contigs, speciesThreshold(refs) ] }
            .combine( SUBTYPE_SKETCH.out.signatures.map { meta, sigs -> [ meta.id, sigs ] }, by: 0 )
            .map { key, meta, contigs, threshold, sigs -> [ meta, contigs, sigs, threshold ] }
    )
    ch_versions = ch_versions.mix(SUBTYPE.out.versions.first())

    // -------------------------------------------------------------------------
    // Match each sample to the reference for its subtype:
    // [ meta, reads, [ matching reference records ] ]
    // -------------------------------------------------------------------------
    ch_species.needs_subtype
        .map { meta, reads, contigs, refs -> [ meta, reads, refs ] }
        .join(
            SUBTYPE.out.subtype
                .splitCsv(header: true)
                .map { meta, row ->
                    def subtype = row['subtype']?.trim()
                    [ meta, subtype && subtype != 'undefined' ? subtype : null ]
                }
        )
        .map { meta, reads, refs, subtype -> [ meta + [ subtype: subtype ], reads, refs ] }
        .mix( ch_species.has_subtype.map { meta, reads, contigs, refs -> [ meta, reads, refs ] } )
        .map { meta, reads, refs ->
            [ meta, reads, refs.findAll { ref -> matchesSubtype(ref, meta.subtype) } ]
        }
        .branch { meta, reads, refs ->
            matched:   refs
            unmatched: true
        }
        .set { ch_matched }

    // -------------------------------------------------------------------------
    // Fail the pipeline, listing every sample without a reference and why
    // -------------------------------------------------------------------------
    ch_species.no_species
        .map { meta, reads, contigs, refs ->
            "${meta.id}: GAMBIT did not return a species-level call"
        }
        .mix(
            ch_species.no_reference.map { meta, reads, contigs, refs ->
                "${meta.id}: no reference for species '${meta.species}'"
            },
            ch_matched.unmatched.map { meta, reads, refs ->
                Utils.hasValue(meta.subtype)
                    ? "${meta.id}: no reference for species '${meta.species}' with subtype '${meta.subtype}'"
                    : "${meta.id}: subtype could not be determined for species '${meta.species}'"
            }
        )
        .map { msg -> "  - ${msg}".toString() }
        .collect()
        .map { lines ->
            error("No matching reference found for ${lines.size()} sample(s):\n" + lines.sort().join('\n'))
        }

    // One item per sample/reference pair: [ meta, reads, ref_meta, reference ]
    ch_matched.matched
        .flatMap { meta, reads, refs ->
            refs.collect { ref -> [ meta, reads, ref.findAll { k, v -> k != 'reference' }, ref.reference ] }
        }
        .set { ch_sample_refs }

    emit:
    samplesheet = ch_matched.matched.map { meta, reads, refs -> [ meta, reads ] }  // species / subtype set
    sample_refs = ch_sample_refs  // [ meta, reads, ref_meta, reference ]
    denovo      = SHOVILL.out.contigs
    species     = GAMBIT_QUERY.out.taxa
    subtype     = SUBTYPE.out.subtype
    signatures  = SUBTYPE_SKETCH.out.signatures  // [ [id: species_key], sig.zip ]
    versions    = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Species from a GAMBIT summary CSV, following subtyper.py: the first row with
// a predicted.name, accepted only when predicted.rank is 'species'.
//
def gambitSpecies(csv) {
    def row = csv.splitCsv(header: true).find { it['predicted.name']?.trim() }
    def isSpecies = row && row['predicted.rank']?.trim()?.toLowerCase() == 'species'
    return isSpecies ? row['predicted.name'].trim() : null
}

//
// Whether a reference record matches a species or subtype name. Names are
// compared after Utils.sanitize(), and an unset name never matches.
//
def matchesSpecies(ref, species) {
    return species ? ref.species.any { Utils.sanitize(it) == Utils.sanitize(species.toString()) } : false
}

def matchesSubtype(ref, subtype) {
    return subtype ? ref.subtype.any { Utils.sanitize(it) == Utils.sanitize(subtype.toString()) } : false
}

//
// Key that groups samples by their species' reference set: the manifest
// entry's name (e.g. 'candidozyma_auris'), so aliases like 'Candida auris'
// share one set of signatures.
//
def speciesKey(refs) {
    return refs[0].name
}

//
// ANI threshold for a species, from the subtype_ani of its reference records,
// as a fraction (e.g. 0.997). subtyper.py takes one threshold per run, so if
// the species' subtypes have different values the strictest (highest) is
// used. Returns '' when none is set, so the script's default applies.
//
def speciesThreshold(refs) {
    def values = refs
        .collect { ref -> ref.subtype_ani }
        .findAll { v -> v != null && v.toString().trim() }
        .collect { v -> v.toString().trim() as BigDecimal }
    return values ? values.max().stripTrailingZeros().toPlainString() : ''
}

//
// Reference record without its path, plus the staged file name, so
// SUBTYPE_SKETCH can tell which assembly belongs to which subtype.
//
def referenceInfo(ref) {
    return ref.findAll { k, v -> k != 'reference' } + [ assembly: ref.reference.name ]
}
