//
// Call consensus peaks with BEDTools and custom scripts, annotate with HOMER, quantify with featureCounts and QC with DESeq2
//

include { HOMER_ANNOTATEPEAKS    } from '../../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS  } from '../../modules/nf-core/subread/featurecounts/main'

include { MACS2_CONSENSUS        } from '../../modules/local/macs2_consensus'
include { COUNT_NORM              } from '../../modules/local/count_normalization' 

workflow BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 {
    take:
    ch_peaks                            // channel: [ val(meta), [ peaks ] ]
    ch_bams                             // channel: [ val(meta), [ bams ] ]
    ch_fasta                            // channel: [ fasta ]
    ch_gtf                              // channel: [ gtf ]
    ch_deseq2_pca_header_multiqc        // channel: [ header_file ]
    ch_deseq2_clustering_header_multiqc // channel: [ header_file ]
    is_narrow_peak                      // boolean: true/false
    skip_peak_annotation                // boolean: true/false
    skip_deseq2_qc                      // boolean: true/false
    
    main:

    ch_versions = Channel.empty()

    // Create channels: [ meta , [ peaks ] ]
    // where meta = [ id : consensus_peaks ]
    ch_peaks
        .collect { it[1] }
        .filter { it.size() > 1 }
        .map { 
            peaks ->
                [ [ id: 'consensus_peaks' ], peaks ]
        }
        .set { ch_consensus_peaks }

    //
    // Generate consensus peaks across samples
    //
    MACS2_CONSENSUS (
        ch_consensus_peaks,
        is_narrow_peak
    )
    ch_versions = ch_versions.mix(MACS2_CONSENSUS.out.versions)

    //
    // Annotate consensus peaks
    //
    ch_homer_annotatepeaks = Channel.empty()
    if (!skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS (
            MACS2_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    }

    // Create channels: [ meta, [ bams ], saf ]
    ch_bams
        .join(ch_peaks)
        .collect { it[1] }
        .filter { it.size() > 1 }
        .map { [ it ] }
        .concat(MACS2_CONSENSUS.out.saf)
        .collect()
        .filter { it.size() == 3 }
        .map {
            bam, meta, saf -> 
                [ meta, bam , saf ]
        }
        .set { ch_bam_saf }

    //
    // Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_bam_saf
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    //
    // Generate QC plots with DESeq2 
    //
    ch_deseq2_qc_pdf           = Channel.empty()
    ch_deseq2_qc_rdata         = Channel.empty()
    ch_deseq2_qc_rds           = Channel.empty()
    ch_deseq2_qc_pca_txt       = Channel.empty()
    ch_deseq2_qc_pca_multiqc   = Channel.empty()
    ch_deseq2_qc_dists_txt     = Channel.empty()
    ch_deseq2_qc_dists_multiqc = Channel.empty()
    ch_deseq2_qc_log           = Channel.empty()
    ch_deseq2_qc_size_factors  = Channel.empty()
    if (!skip_deseq2_qc) {
        COUNT_NORM (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = COUNT_NORM.out.pdf
        ch_deseq2_qc_rdata         = COUNT_NORM.out.rdata
        ch_deseq2_qc_rds           = COUNT_NORM.out.rds
        ch_deseq2_qc_pca_txt       = COUNT_NORM.out.pca_txt
        ch_deseq2_qc_pca_multiqc   = COUNT_NORM.out.pca_multiqc
        ch_deseq2_qc_dists_txt     = COUNT_NORM.out.dists_txt
        ch_deseq2_qc_dists_multiqc = COUNT_NORM.out.dists_multiqc
        ch_deseq2_qc_log           = COUNT_NORM.out.log
        ch_deseq2_qc_size_factors  = COUNT_NORM.out.size_factors

        ch_deseq2_qc_noamlization       = COUNT_NORM.out.noamlization
        ch_deseq2_qc_noamlization_txt   = COUNT_NORM.out.noamlization_txt


        ch_versions = ch_versions.mix(COUNT_NORM.out.versions)
    }

    emit:
    consensus_bed           = MACS2_CONSENSUS.out.bed           // channel: [ bed ]
    consensus_saf           = MACS2_CONSENSUS.out.saf           // channel: [ saf ]
    consensus_pdf           = MACS2_CONSENSUS.out.pdf           // channel: [ pdf ]
    consensus_boolean_txt   = MACS2_CONSENSUS.out.boolean_txt   // channel: [ txt ]
    consensus_intersect_txt = MACS2_CONSENSUS.out.intersect_txt // channel: [ txt ]

    homer_annotatepeaks     = ch_homer_annotatepeaks            // channel: [ txt ]

    featurecounts_txt       = SUBREAD_FEATURECOUNTS.out.counts  // channel: [ txt ]
    featurecounts_summary   = SUBREAD_FEATURECOUNTS.out.summary // channel: [ txt ]

    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_pca_multiqc   = ch_deseq2_qc_pca_multiqc          // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_dists_multiqc = ch_deseq2_qc_dists_multiqc        // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]
    deseq2_qc_noamlization  = ch_deseq2_qc_noamlization
    deseq2_qc_noamlization_txt  = ch_deseq2_qc_noamlization_txt

    versions                = ch_versions                       // channel: [ versions.yml ]
}
