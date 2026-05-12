// modules/fastqc.nf
// Raw read quality control with FastQC

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    conda     'envs/rnaseq.yml'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.html'), emit: html
    tuple val(meta), path('*.zip'),  emit: zip
    path '*.version.txt',            emit: version

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${args} \\
        ${reads}

    fastqc --version | sed 's/FastQC //' > fastqc.version.txt
    """

    stub:
    """
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip fastqc.version.txt
    """
}
