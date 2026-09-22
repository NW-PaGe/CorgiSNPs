/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_corgisnps_pipeline'

// Subworkflows
include { PREPARE  } from '../subworkflows/local/prepare'
include { CLASSIFY } from '../subworkflows/local/classify'
include { VARIANTS } from '../subworkflows/local/variants'
include { AMR      } from '../subworkflows/local/amr'
include { PHYLO    } from '../subworkflows/local/phylo'

// Report / summary modules
include { SUMMARYLINE    } from '../modules/local/report/main'
include { ADD_AMR        } from '../modules/local/report/main'
include { REPORT_SPECIES } from '../modules/local/report/main'
include { REPORT_ALL     } from '../modules/local/report/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow CORGISNPS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // =========================================================================
    // PREPARE
    // =========================================================================
    PREPARE(ch_samplesheet)

    ch_versions      = ch_versions.mix(PREPARE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE.out.multiqc_files)
    ch_refs          = PREPARE.out.refs
    ch_samplesheet   = PREPARE.out.samplesheet
    ch_read_stats    = PREPARE.out.read_stats

    // Placeholder channels ([meta, []]) for optional downstream outputs
    ch_blank     = ch_samplesheet.map { meta, reads -> [meta, []] }
    ch_denovo    = ch_blank
    ch_species   = ch_blank
    ch_subtype   = ch_blank
    ch_aln_stats = ch_blank
    ch_tree      = ch_blank
    ch_dist      = ch_blank

    // =========================================================================
    // CLASSIFY (optional) - only samples missing species or subtype
    // =========================================================================
    if (params.classify) {
        ch_to_classify   = ch_samplesheet.filter { meta, reads -> !(meta.species && meta.subtype) }
        ch_pre_classified = ch_samplesheet.filter { meta, reads ->   meta.species && meta.subtype  }

        CLASSIFY(ch_to_classify, ch_refs)

        ch_versions    = ch_versions.mix(CLASSIFY.out.versions)
        ch_samplesheet = CLASSIFY.out.samplesheet.concat(ch_pre_classified)
        ch_denovo      = CLASSIFY.out.denovo
        ch_species     = CLASSIFY.out.species
        ch_subtype     = CLASSIFY.out.subtype
    }

    // Sanitize species/subtype strings
    ch_samplesheet = ch_samplesheet.map { meta, reads ->
        [ meta + [species: Utils.sanitize(meta.species), subtype: Utils.sanitize(meta.subtype)], reads ]
    }

    // =========================================================================
    // SUMMARYLINE - per-sample summary and auto QC
    // Joins keep samples lacking some inputs (remainder: true); missing
    // entries become [] so every tuple has the same shape.
    // =========================================================================
    ch_samples = ch_samplesheet
        .map  { meta, reads -> [meta.id, meta] }
        .join (ch_read_stats.map { meta, file -> [meta.id, file] }, remainder: true)
        .join (ch_denovo.map     { meta, file -> [meta.id, file] }, remainder: true)
        .join (ch_species.map    { meta, file -> [meta.id, file] }, remainder: true)
        .join (ch_subtype.map    { meta, file -> [meta.id, file] }, remainder: true)
        .map  { row -> row.drop(1).collect { it ?: [] } } // drop join key, replace nulls with []

    SUMMARYLINE(
        ch_samples,
        file(params.input),
        file(params.ncbi_stats)
    )
    ch_versions = ch_versions.mix(SUMMARYLINE.out.versions.first())

    // QC status per sample (--ignore_qc passes everything)
    ch_auto_qc = SUMMARYLINE.out.summary
        .splitCsv(header: true)
        .map    { meta, data -> [meta, params.ignore_qc || data.qc_status == 'PASS'] }
        .branch { meta, status ->
            pass    : status
            not_pass: !status
        }

    ch_qc_pass_meta = ch_auto_qc.pass.map     { meta, status -> [meta] }
    ch_qc_fail_meta = ch_auto_qc.not_pass.map { meta, status -> [meta] }

    ch_samplesheet_pass = ch_samplesheet.join(ch_qc_pass_meta)
    ch_summary_pass     = SUMMARYLINE.out.summary.join(ch_qc_pass_meta)
    ch_summary_fail     = SUMMARYLINE.out.summary.join(ch_qc_fail_meta)

    // =========================================================================
    // VARIANTS / AMR / PHYLO (optional)
    // =========================================================================
    if (params.variants) {
        VARIANTS(ch_samplesheet_pass, ch_refs, true)

        ch_versions         = ch_versions.mix(VARIANTS.out.versions)
        ch_samplesheet_pass = VARIANTS.out.samplesheet
        ch_bam              = VARIANTS.out.bam
        ch_vcf              = VARIANTS.out.vcf
        ch_aln              = VARIANTS.out.aln

        if (params.amr) {
            AMR(ch_samplesheet_pass, ch_bam, ch_vcf, ch_refs)
            ch_versions = ch_versions.mix(AMR.out.versions)

            // Join summary rows with AMR output by sample ID (meta.id), not full meta
            ch_summary_pass_amr = ch_summary_pass
                .map  { meta, summaryline -> [meta.id, meta, summaryline] }
                .join (AMR.out.summary.map { meta, amr_summary -> [meta.id, amr_summary] }, remainder: true)
                .map  { id, meta, summaryline, amr_summary -> [meta, summaryline, amr_summary] }
                .branch { meta, summaryline, amr_summary ->
                    pass    : summaryline && amr_summary
                    not_pass: summaryline && !amr_summary
                }

            ADD_AMR(ch_summary_pass_amr.pass)
            ch_versions = ch_versions.mix(ADD_AMR.out.versions)

            ch_summary_pass = ADD_AMR.out.summary.concat(
                ch_summary_pass_amr.not_pass.map { meta, summaryline, amr_summary -> [meta, summaryline] }
            )
        }

        if (params.phylo) {
            PHYLO(ch_aln, ch_samplesheet_pass)

            ch_versions  = ch_versions.mix(PHYLO.out.versions)
            ch_aln_stats = PHYLO.out.aln_stats
            ch_tree      = PHYLO.out.tree
            ch_dist      = PHYLO.out.dist
        }
    }

    // =========================================================================
    // REPORTS
    // =========================================================================
    REPORT_ALL(
        ch_summary_pass
            .concat(ch_summary_fail)
            .map { meta, summaryline -> summaryline }
            .collect()
    )

    if (params.phylo) {
        REPORT_SPECIES(
            ch_aln_stats
                .join(ch_dist, by: 0)
                .join(ch_tree, by: 0, remainder: true)
                .map { meta, aln_stats, dist, tree -> [meta, aln_stats, tree ?: [], dist] }
                .combine(REPORT_ALL.out.summary),
            file(params.microreact_template)
        )
    }

    // =========================================================================
    // Collate and save software versions
    // =========================================================================
    ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name    : 'CorgiSNPs_software_mqc_versions.yml',
            sort    : true,
            newLine : true
        )

    // =========================================================================
    // MultiQC
    // =========================================================================
    ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo,   checkIfExists: true) : Channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_methods_description_file = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_methods_description_file))

    ch_multiqc_files = ch_multiqc_files
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // path: MultiQC HTML report
    versions       = ch_versions                 // channel: versions.yml files from all stages
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
