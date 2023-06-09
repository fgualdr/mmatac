process BAM_SIEVE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.shifted.bam"), emit: bam_shift
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? task.memory.toGiga() : false
    def sort_mem = avail_mem ? "-m ${
                                        def mem_gb = Math.round((avail_mem - 4) / task.cpus)
                                        mem_gb == 0 ? '2G' : "${mem_gb}G"
                                    }" : '-m 2G'
    """
    # 1) Shift the 5' end +4 bp for positive and -5bp for negatively mapped reads.
    samtools view -h $bam | awk 'BEGIN{OFS="\t"} $1~/^@/{print;next} {if($9>0){$4=$4+4;print}else{$4=$4-5;print}}' | samtools view -bS - > ${prefix}.shifted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
        deeptools: $(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """


}