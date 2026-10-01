//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { POLYCORE } from '../../../modules/local/polycore/main'
include { IQTREE   } from '../../../modules/local/iqtree/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PHYLO {

    take:
    // ch_aln : [ val(meta), path(aln) ]  // consensus alignments (label typo kept)
    ch_aln
    // ch_samplesheet: [ val(meta), path(reads) ]
    ch_samplesheet

    main:

    // Collector for version files
    ch_versions = Channel.empty()

    // -------------------------------------------------------------------------
    // Join alignments with meta, group by (species, subtype), and assemble a
    // sample set consisting of current alignments plus any cached consensus
    // sequences from params.db not already present in the current set.
    //
    // Groups are keyed on species/subtype/reference/ploidy only, so samples
    // with different settings for the same reference still end up in one
    // tree. The group's phylo settings are those its samples agree on; if
    // they disagree, the run-level params are used, with a warning.
    // -------------------------------------------------------------------------
    ch_sp_sb = ch_aln
        .map{ meta, aln ->
            def key = [species: meta.species, subtype: meta.subtype, reference: meta.reference, ploidy: meta.ploidy]
            return [key, meta.settings ?: [:], aln]
        }
        .groupTuple(by: 0)
        .map { key, settings, alns ->
            def meta = key + [ settings: groupSettings(key, settings) ]
            // Current sample names (filenames) from the alignment paths
            def current = alns.collect { it.getName() }

            // Collect cached consensus sequences for this (sp/sb) from params.db
            def filtered_cache = []
            def db_path        = file(params.db)
            if (db_path.exists()){
                def sp_sb_db    = db_path.resolve(meta.species).resolve(meta.subtype)
                def sp_db_cache = sp_sb_db.exists() ? sp_sb_db.list().collect { sp_sb_db.resolve(it) } : []
                filtered_cache  = sp_db_cache.findAll { !current.contains(it.getName()) }
            }
            

            // Combine current alignments with cache, keeping ref as a single path
            def all_alns = alns + filtered_cache
            [ meta, meta.reference, all_alns ]
        }

    // -------------------------------------------------------------------------
    // POLYCORE: derive SNP alignment, distance matrix, and per-sample stats
    // Input: [ species, subtype, samples, reference ]
    // -------------------------------------------------------------------------
    POLYCORE(
        ch_sp_sb
    )
    ch_versions = ch_versions.mix(POLYCORE.out.versions.first())

    // -------------------------------------------------------------------------
    // IQTREE: build trees only when there are >2 samples for the clade
    // Joins: SNP alignment + constant-site file + sample count
    // -------------------------------------------------------------------------
    IQTREE(
        POLYCORE.out.snps
            .join(POLYCORE.out.fconst)
            .join(POLYCORE.out.uniq_seq)
            .filter { meta, aln, const_sites, uniq_seq ->
                if( (uniq_seq?.isInteger() ? uniq_seq.toInteger() : 0) > 2 ){
                    return true
                } else {
                    log.warn "${meta.species} ${meta.subtype} contains <3 unique sequences - no tree will be created!"
                    return false
                }
            }
    )
    ch_versions = ch_versions.mix(IQTREE.out.versions.first())

    emit:
    aln_stats   = POLYCORE.out.csv
    tree        = IQTREE.out.tree
    dist        = POLYCORE.out.dist_wide
    versions    = ch_versions
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Phylogenetic settings for a species/subtype group, from its samples'
// meta.settings. Only phylo settings are kept (variant calling settings are
// per sample), so the group meta is the same however the samples arrive.
// If the samples disagree, the run-level params are used.
//
def groupSettings(key, List settings) {
    def phyloKeys = ReferenceManifest.PHYLO_SETTINGS.keySet()
    def distinct  = settings.collect { s -> s.subMap(phyloKeys) }.unique()
    if( distinct.size() > 1 )
        log.warn "${key.species} ${key.subtype}: samples have different phylogenetic settings; using the run-level parameters for this group"
    def chosen = distinct.size() == 1 ? distinct[0] : [:]
    return ReferenceManifest.resolveSettings(chosen, params).subMap(phyloKeys)
}
