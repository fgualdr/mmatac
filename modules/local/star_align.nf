process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "bioconda::star=2.7.10a bioconda::samtools=1.16.1 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'quay.io/biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path "versions.yml"                       , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    
    def out_sam_type = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'    // this is the default
    
    def filtermultimapnmax = (args.contains('--outFilterMultimapNmax')) ? '' : "--outFilterMultimapNmax $params.outfiltermultimapnmax"   
    def outsammultinmax = (args.contains('--outSAMmultNmax')) ? '' : "--outSAMmultNmax $params.outsammultnmax"
    def anchormultimapnmax = (args.contains('--winAnchorMultimapNmax')) ? '' : "--winAnchorMultimapNmax $params.winanchormultimapnmax"

    def outmultiorder = (args.contains('--outMultimapperOrder')) ? '' : '--outMultimapperOrder Random'
    def endtoendprotrude = (args.contains('--alignEndsProtrude')) ? '' : '--alignEndsProtrude 0 ConcordantPair' 
    def endtoendtype = (args.contains('--alignEndsType')) ? '' : '--alignEndsType EndToEnd' 
    def mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --outFilterScoreMinOverLread 0.3 \\
        --outFilterMatchNminOverLread 0.3 \\
        --outSAMattrIHstart 0 \\
        $filtermultimapnmax \\
        $outsammultinmax \\
        $anchormultimapnmax \\
        $outmultiorder \\
        $endtoendprotrude \\
        $endtoendtype \\
        $out_sam_type \\
        $seq_center \\
        $args

    $mv_unsorted_bam
    
    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}
