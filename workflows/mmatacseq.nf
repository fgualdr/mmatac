/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners : [ 'star' ]
]
 
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMMAtacseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta,
    params.gtf, params.gff, params.gene_bed, params.tss_bed,
    params.star_index,
    params.blacklist,
    params.bamtools_filter_pe_config, params.bamtools_filter_se_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


// JSON files required by BAMTools for alignment filtering
ch_bamtools_filter_se_config = file(params.bamtools_filter_se_config, checkIfExists: true)
ch_bamtools_filter_pe_config = file(params.bamtools_filter_pe_config, checkIfExists: true)

// Header files for MultiQC
ch_multiqc_merged_library_peak_count_header        = file("$projectDir/assets/multiqc/merged_library_peak_count_header.txt", checkIfExists: true)
ch_multiqc_merged_library_frip_score_header        = file("$projectDir/assets/multiqc/merged_library_frip_score_header.txt", checkIfExists: true)
ch_multiqc_merged_library_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_library_peak_annotation_header.txt", checkIfExists: true)
ch_multiqc_merged_library_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_library_deseq2_pca_header.txt", checkIfExists: true)
ch_multiqc_merged_library_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_library_deseq2_clustering_header.txt", checkIfExists: true)

ch_multiqc_merged_replicate_peak_count_header        = file("$projectDir/assets/multiqc/merged_replicate_peak_count_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_frip_score_header        = file("$projectDir/assets/multiqc/merged_replicate_frip_score_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_peak_annotation_header   = file("$projectDir/assets/multiqc/merged_replicate_peak_annotation_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_deseq2_pca_header        = file("$projectDir/assets/multiqc/merged_replicate_deseq2_pca_header.txt", checkIfExists: true)
ch_multiqc_merged_replicate_deseq2_clustering_header = file("$projectDir/assets/multiqc/merged_replicate_deseq2_clustering_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IGV     } from '../modules/local/igv'
include { MULTIQC } from '../modules/local/multiqc'
include { MERGED_LIBRARY_DEEPTOOL_BIGWIG   } from '../modules/local/merged_library_deep_bigwig'
include { MERGED_LIBRARY_DEEPTOOLS_BIGWIG_NORM   } from '../modules/local/merged_library_deep_bigwig_norm'
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR     } from '../subworkflows/local/align_star'
include { BIGWIG_PLOT_DEEPTOOLS as MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS       } from '../subworkflows/local/bigwig_plot_deeptools'
include { BAM_FILTER_BAMTOOLS as MERGED_LIBRARY_FILTER_BAM                    } from '../subworkflows/local/bam_filter_bamtools'



include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS2_HOMER as MERGED_LIBRARY_CALL_ANNOTATE_PEAKS   } from '../subworkflows/local/bam_peaks_call_qc_annotate_macs2_homer.nf'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 as MERGED_LIBRARY_CONSENSUS_PEAKS   } from '../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PICARD_COLLECTMULTIPLEMETRICS as MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP as MERGED_LIBRARY_PRESEQ_LCEXTRAP                             } from '../modules/nf-core/preseq/lcextrap/main'
include { DEEPTOOLS_PLOTFINGERPRINT as MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT         } from '../modules/nf-core/deeptools/plotfingerprint/main'
include { ATAQV_ATAQV as MERGED_LIBRARY_ATAQV_ATAQV                                     } from '../modules/nf-core/ataqv/ataqv/main'
include { ATAQV_MKARV as MERGED_LIBRARY_ATAQV_MKARV                                     } from '../modules/nf-core/ataqv/mkarv/main'

include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_LIBRARY   } from '../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_MERGESAMFILES as PICARD_MERGESAMFILES_REPLICATE } from '../modules/nf-core/picard/mergesamfiles/main'

include { MACS2_CALLPEAK   as    MACS2_CALLPEAK_MERGED     } from '../modules/nf-core/macs2/callpeak/main'

//
// SUBWORKFLOW: 
//
include { FASTQ_FASTQC_UMITOOLS_FASTP } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'

include { BAM_MARKDUPLICATES_PICARD as MERGED_LIBRARY_MARKDUPLICATES_PICARD   } from '../subworkflows/nf-core/bam_markduplicates_picard/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MMATACSEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.aligner
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            INPUT_CHECK.out.reads,
            params.skip_fastqc || params.skip_qc,
            false,
            false,
            params.skip_trimming,
            0
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            INPUT_CHECK.out.reads,
            params.skip_fastqc || params.skip_qc,
            false,
            false,
            0,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    
    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    ch_star_multiqc = Channel.empty()
    if (params.aligner == 'star') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.fasta,
            params.seq_center ?: ''
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)
    }

    // Create channels: [ meta, [bam] ]
    ch_genome_bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id - ~/_T\d+$/
                [ meta_clone, bam ]
        }
        .groupTuple(by: [0])
        .map {
            meta, bam ->
                [ meta, bam.flatten() ]
        }
        .set { ch_sort_bam }

    //
    // MODULE: Merge resequenced BAM files
    //
    PICARD_MERGESAMFILES_LIBRARY (
        ch_sort_bam
    )
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES_LIBRARY.out.versions.first())

    //
    // SUBWORKFLOW: Mark duplicates in BAM files
    //
    MERGED_LIBRARY_MARKDUPLICATES_PICARD (
        PICARD_MERGESAMFILES_LIBRARY.out.bam,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.versions)

    //
    // SUBWORKFLOW: Filter BAM file
    //
    MERGED_LIBRARY_FILTER_BAM (
        MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bam.join(MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.bai, by: [0]),
        PREPARE_GENOME.out.filtered_bed.first(),
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.chrom_sizes,
        ch_bamtools_filter_se_config,
        ch_bamtools_filter_pe_config
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_FILTER_BAM.out.versions)

    //
    // MODULE: Picard post alignment QC
    //
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS (
            MERGED_LIBRARY_FILTER_BAM.out.bam,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
        )
        ch_picardcollectmultiplemetrics_multiqc = MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(MERGED_LIBRARY_PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }


    // Create channels: [ meta, [bam], [bai] ]
    MERGED_LIBRARY_FILTER_BAM
        .out
        .bam
        .join(MERGED_LIBRARY_FILTER_BAM.out.bai, by: [0])
        .set { ch_bam_bai }

    //
    // SUBWORKFLOW: BigWig coverage tracks CPM with deeptools:
    //
    // MERGED_LIBRARY_DEEPTOOL_BIGWIG (
    //     ch_bam_bai
    // )
    // ch_versions = ch_versions.mix(MERGED_LIBRARY_DEEPTOOL_BIGWIG.out.versions)

    //
    // MODULE: deepTools plotFingerprint QC
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint) {
        MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT (
            ch_bam_bai
        )
        ch_deeptoolsplotfingerprint_multiqc = MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    // Create channels: [ meta, bam, ([] for control_bam) ]
    ch_bam_bai
        .map {
            meta, bam, bai ->
                [ meta , bam, [] ]
        }
        .set { ch_bam_library }

    // Call Peaks on merged BAMS:
    ch_bam_bai
        .map {
            meta, bam, bai ->
            def exp = meta.experiment
            def new_meta = meta.clone()
            new_meta.id =  meta.experiment
            [new_meta, bam, bai]
        }
        .groupTuple(by: 0)
        .map {
            meta, bam, bai ->
                [ meta , bam, [] ]
        }
        .set { ch_bam_exp }

    //
    // S3norm pipe yes/no
    // 

    // 
    // If we don't use the S3norm pipe:
    MACS2_CALLPEAK_MERGED(
        ch_bam_exp,
        PREPARE_GENOME.out.macs_gsize
    )
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_MERGED.out.versions.first())
    //
    // SUBWORKFLOW: Call peaks with MACS2, annotate with HOMER and perform downstream QC
    // 
    MERGED_LIBRARY_CALL_ANNOTATE_PEAKS (
        ch_bam_library,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.macs_gsize,
        ".mLb.clN_peaks.annotatePeaks.txt",
        ch_multiqc_merged_library_peak_count_header,
        ch_multiqc_merged_library_frip_score_header,
        ch_multiqc_merged_library_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.versions)
    //
    // SUBWORKFLOW: Consensus peaks analysis
    //
    ch_macs2_consensus_library_bed       = Channel.empty()
    ch_featurecounts_library_multiqc     = Channel.empty()
    ch_deseq2_pca_library_multiqc        = Channel.empty()
    ch_deseq2_clustering_library_multiqc = Channel.empty()

    MERGED_LIBRARY_CONSENSUS_PEAKS (
        MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peaks,
        ch_bam_library,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.gtf,
        ch_multiqc_merged_library_deseq2_pca_header,
        ch_multiqc_merged_library_deseq2_clustering_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_deseq2_qc
    )
    ch_macs2_consensus_library_bed       = MERGED_LIBRARY_CONSENSUS_PEAKS.out.consensus_bed
    ch_featurecounts_library_multiqc     = MERGED_LIBRARY_CONSENSUS_PEAKS.out.featurecounts_summary
    ch_deseq2_pca_library_multiqc        = MERGED_LIBRARY_CONSENSUS_PEAKS.out.deseq2_qc_pca_multiqc
    ch_deseq2_clustering_library_multiqc = MERGED_LIBRARY_CONSENSUS_PEAKS.out.deseq2_qc_dists_multiqc
    ch_versions = ch_versions.mix(MERGED_LIBRARY_CONSENSUS_PEAKS.out.versions)

    MERGED_LIBRARY_CONSENSUS_PEAKS
        .out
        .deseq2_qc_noamlization_txt
        .splitCsv ( header:true, sep:'\t' )
        .map { row -> 
            def id = row.Sample_ID
            def value = row.scaling
            [ id, value ]
        }
        .set { ch_size_factors }
    ch_size_factors.view()

    //
    // Make the normalized BigWigs
    // BigWigs
    //
    ch_bam_bai
        .combine(ch_size_factors)
        .map { 
            meta1, bam1, bai1, id2, scaling2 ->
                meta1.id == id2 ? [ meta1, bam1, bai1 ,scaling2] : null
        }
        .set { ch_bam_bai_scale } 
    ch_bam_bai_scale.view()
    //
    // MODULE: BedGraph coverage tracks.
    //
    MERGED_LIBRARY_DEEPTOOLS_BIGWIG_NORM (
         ch_bam_bai_scale // join with the created channel of scalings - Use deeptool insrtead make directly wiggles
    )
    ch_versions = ch_versions.mix(MERGED_LIBRARY_DEEPTOOLS_BIGWIG_NORM.out.versions.first())

    // 
    // Use S3norm
    //

    //
    // SUBWORKFLOW: Plot coverage across annotation with deepTools
    //
    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    if (!params.skip_plot_profile) {
        MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS (
            MERGED_LIBRARY_DEEPTOOLS_BIGWIG_NORM.out.bigwig,
            PREPARE_GENOME.out.gene_bed,
            PREPARE_GENOME.out.tss_bed
        )
        ch_deeptoolsplotprofile_multiqc = MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS.out.plotprofile_table
        ch_versions = ch_versions.mix(MERGED_LIBRARY_BIGWIG_PLOT_DEEPTOOLS.out.versions)
    }


    // Create channels: [ meta, bam, bai, peak_file ]
    MERGED_LIBRARY_FILTER_BAM
        .out
        .bam
        .join(MERGED_LIBRARY_FILTER_BAM.out.bai, by: [0])
        .join(MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peaks, by: [0])
        .set { ch_bam_peaks }

    //
    // MODULE: ATAQV QC
    //
    if (!params.skip_ataqv) {
        MERGED_LIBRARY_ATAQV_ATAQV (
            ch_bam_peaks,
            'NA',
            params.mito_name,
            PREPARE_GENOME.out.tss_bed,
            [],
            PREPARE_GENOME.out.autosomes
        )
        ch_versions = ch_versions.mix(MERGED_LIBRARY_ATAQV_ATAQV.out.versions.first())

        MERGED_LIBRARY_ATAQV_MKARV (
            MERGED_LIBRARY_ATAQV_ATAQV.out.json.collect{it[1]}
        )
        ch_versions = ch_versions.mix(MERGED_LIBRARY_ATAQV_MKARV.out.versions)
    }

    //
    // Merged replicate analysis
    //
    ch_markduplicates_replicate_stats                   = Channel.empty()
    ch_markduplicates_replicate_flagstat                = Channel.empty()
    ch_markduplicates_replicate_idxstats                = Channel.empty()
    ch_markduplicates_replicate_metrics                 = Channel.empty()
    ch_ucsc_bedgraphtobigwig_replicate_bigwig           = Channel.empty()
    ch_macs2_replicate_peaks                            = Channel.empty()
    ch_macs2_frip_replicate_multiqc                     = Channel.empty()
    ch_macs2_peak_count_replicate_multiqc               = Channel.empty()
    ch_macs2_plot_homer_annotatepeaks_replicate_multiqc = Channel.empty()
    ch_macs2_consensus_replicate_bed                    = Channel.empty()
    ch_featurecounts_replicate_multiqc                  = Channel.empty()
    ch_deseq2_pca_replicate_multiqc                     = Channel.empty()
    ch_deseq2_clustering_replicate_multiqc              = Channel.empty()
    ch_preseq_multiqc                                   = Channel.empty()

    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {

        workflow_summary    = WorkflowMMAtacseq.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)
        methods_description    = WorkflowMMAtacseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

            ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]),
            ch_trim_log_multiqc.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.stats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.flagstat.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.idxstats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_MARKDUPLICATES_PICARD.out.metrics.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_FILTER_BAM.out.stats.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_FILTER_BAM.out.flagstat.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_FILTER_BAM.out.idxstats.collect{it[1]}.ifEmpty([]),
            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),

            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),

            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.frip_multiqc.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.peak_count_multiqc.collect{it[1]}.ifEmpty([]),
            MERGED_LIBRARY_CALL_ANNOTATE_PEAKS.out.plot_homer_annotatepeaks_tsv.collect().ifEmpty([]),
            ch_featurecounts_library_multiqc.collect{it[1]}.ifEmpty([]),

            ch_markduplicates_replicate_stats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_flagstat.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_replicate_metrics.collect{it[1]}.ifEmpty([]),

            ch_macs2_frip_replicate_multiqc.collect{it[1]}.ifEmpty([]),
            ch_macs2_peak_count_replicate_multiqc.collect{it[1]}.ifEmpty([]),
            ch_macs2_plot_homer_annotatepeaks_replicate_multiqc.collect().ifEmpty([]),
            ch_featurecounts_replicate_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deseq2_pca_library_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_library_multiqc.collect().ifEmpty([]),
            ch_deseq2_pca_replicate_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_replicate_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }

    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
