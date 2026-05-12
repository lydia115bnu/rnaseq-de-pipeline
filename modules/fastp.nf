// modules/fastp.nf
// Adapter trimming and quality filtering with fastp

process FASTP {
    tag "${meta.id}"
    label 'process_medium'

    conda     'envs/rnaseq.yml'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastp.fastq.gz'), emit: trimmed
    tuple val(meta), path('*.json'),           emit: json
    tuple val(meta), path('*.html'),           emit: html
    path '*.version.txt',                      emit: version

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (r1, r2) = reads

    // Paired-end mode
    """
    fastp \\
        --in1 ${r1} \\
        --in2 ${r2} \\
        --out1 ${prefix}_R1.fastp.fastq.gz \\
        --out2 ${prefix}_R2.fastp.fastq.gz \\
        --json ${prefix}.fastp.json \\
        --html ${prefix}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --correction \\
        --qualified_quality_phred 20 \\
        --length_required 36 \\
        ${args}

    fastp --version 2>&1 | head -1 | sed 's/fastp //' > fastp.version.txt
    """

    stub:
    """
    touch ${meta.id}_R1.fastp.fastq.gz ${meta.id}_R2.fastp.fastq.gz
    touch ${meta.id}.fastp.json ${meta.id}.fastp.html fastp.version.txt
    """
}
