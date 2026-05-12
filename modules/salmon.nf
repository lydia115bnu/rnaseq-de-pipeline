// modules/salmon.nf
// Quasi-mapping quantification with Salmon
// Uses pre-built decoy-aware transcriptome index

process SALMON_INDEX {
    tag "building_index"
    label 'process_high'

    conda     'envs/rnaseq.yml'
    container 'quay.io/biocontainers/salmon:1.10.3--h3dcb392_0'

    input:
    path transcriptome   // e.g. Homo_sapiens.GRCh38.cdna.all.fa.gz
    path genome          // full genome FASTA (used for decoy)

    output:
    path 'salmon_index', emit: index
    path '*.version.txt', emit: version

    script:
    """
    # Build decoy list from genome chromosome names
    grep "^>" ${genome} | cut -d " " -f 1 | sed 's/>//' > decoys.txt

    # Concatenate transcriptome + genome (decoy-aware)
    cat ${transcriptome} ${genome} > gentrome.fa.gz

    salmon index \\
        --threads ${task.cpus} \\
        --transcripts gentrome.fa.gz \\
        --decoys decoys.txt \\
        --index salmon_index \\
        --gencode

    salmon --version 2>&1 | sed 's/salmon //' > salmon.version.txt
    """
}

process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_medium'

    conda     'envs/rnaseq.yml'
    container 'quay.io/biocontainers/salmon:1.10.3--h3dcb392_0'

    publishDir "${params.outdir}/salmon/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"),      emit: quants
    tuple val(meta), path("${meta.id}/logs"), emit: log
    path '*.version.txt',                     emit: version

    script:
    def args   = task.ext.args ?: ''
    def (r1, r2) = reads
    """
    salmon quant \\
        --index ${index} \\
        --libType A \\
        --mates1 ${r1} \\
        --mates2 ${r2} \\
        --validateMappings \\
        --gcBias \\
        --seqBias \\
        --numBootstraps 100 \\
        --threads ${task.cpus} \\
        --output ${meta.id} \\
        ${args}

    salmon --version 2>&1 | sed 's/salmon //' > salmon.version.txt
    """

    stub:
    """
    mkdir -p ${meta.id}/logs
    touch ${meta.id}/quant.sf ${meta.id}/logs/salmon_quant.log salmon.version.txt
    """
}
