//
// Subworkflow with functionality specific to the DOH-JDJ0303/mycosnp pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTERQDUMP  } from '../../../modules/local/fasterq-dump/main'
include { SEQTK_SAMPLE } from '../../../modules/local/seqtk/sample/main'
include { FASTQC       } from '../../../modules/nf-core/fastqc/main'
include { FASTP        } from '../../../modules/nf-core/fastp/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPARE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow PREPARE {

    take:
    // Channel: [ meta: meta, reads: reads ] records; meta may include species, subtype, reference, sra
    ch_samplesheet

    main:

    // Collectors
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // -------------------------------------------------------------------------
    // Load the reference manifest (validated unless --validate_refs false).
    // A manifest that can't be read at all is always reported.
    // -------------------------------------------------------------------------
    def ref_db = ReferenceManifest.load(params.reference_db, params.validate_refs as boolean)

    if( ref_db.errors )
        error "Invalid reference directory '${params.reference_db}':\n" + ref_db.errors.collect { "  - ${it}" }.join('\n')

    ch_refs = channel.fromList(ref_db.records)

    // -------------------------------------------------------------------------
    // MODULE: Download reads from SRA for rows with an SRA accession
    // Input reshaped to [ meta, sra ]
    // -------------------------------------------------------------------------
    FASTERQDUMP(
        ch_samplesheet
            .filter{ it.meta.sra }
            .map { it -> [ it.meta, it.meta.sra ] }
    )
    ch_versions = ch_versions.mix(FASTERQDUMP.out.versions)

    // -------------------------------------------------------------------------
    // Merge SRA-derived reads back with pass-through non-SRA reads.
    // Reads are always a list (a single file comes out of a process as a lone
    // path), meta.single_end is set from the number of files for SRA reads,
    // and meta.sra is removed.
    // Output: [ meta: meta, reads: [ paths ] ]
    // -------------------------------------------------------------------------
    FASTERQDUMP
        .out
        .reads
        .map{ meta, reads ->
            def files = reads instanceof List ? reads : [ reads ]
            [ meta: meta + [ single_end: files.size() == 1 ], reads: files ]
        }
        .concat(
            ch_samplesheet.filter{ ! it.meta.sra }
        )
        .map{ it ->
            def files = it.reads instanceof List ? it.reads : [ it.reads ]
            [ meta: it.meta.findAll{ k, v -> k != 'sra' }, reads: files ]
        }
        .set { ch_samplesheet }

    // -------------------------------------------------------------------------
    // MODULE: Downsample reads with seqtk sample (if --max_reads provided)
    // -------------------------------------------------------------------------
    if (params.max_reads) {

        // Compute total read count (approx) by counting first R1 and doubling (paired)
        ch_samplesheet
            .map { it -> [it, it.reads[0].countFastq() * (it.meta.single_end ? 1 : 2) ] }
            .branch { it, n ->
                ok  : n <= params.max_reads
                high: n >  params.max_reads
            }
            .set { ch_reads_count }

        // For high-coverage samples, sample each mate independently
        SEQTK_SAMPLE(
            ch_reads_count
                .high
                .map{ it, n -> [it.meta, it.reads, params.max_reads] }
                .transpose()
        )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)

        // Re-assemble paired reads (sorted by name so R1 comes before R2,
        // since groupTuple collects them in completion order) and merge with
        // the samples that were under the limit
        SEQTK_SAMPLE
            .out
            .read
            .groupTuple(by: 0)
            .map{ meta, reads -> [ meta: meta, reads: reads.sort { r -> r.name } ] }
            .concat( ch_reads_count.ok.map { it, n -> it } )
            .set { ch_samplesheet }
    }

    ch_samplesheet = ch_samplesheet.map{ [it.meta, it.reads] }

    // -------------------------------------------------------------------------
    // MODULE: FastQC (adds zips to MultiQC input and versions to collector)
    // -------------------------------------------------------------------------
    FASTQC(
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect { it[1] } )
    ch_versions      = ch_versions     .mix( FASTQC.out.versions.first() )

    // -------------------------------------------------------------------------
    // MODULE: Fastp (trimming/filters + JSON stats). Replace reads with trimmed.
    // -------------------------------------------------------------------------
    FASTP(
        ch_samplesheet,
        [],
        false,
        false,
        false
    )
    ch_versions    = ch_versions.mix(FASTP.out.versions.first())
    ch_samplesheet = FASTP.out.reads

    emit:
    refs          = ch_refs
    samplesheet   = ch_samplesheet
    read_stats    = FASTP.out.json
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
