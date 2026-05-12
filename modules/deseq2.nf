// modules/deseq2.nf
// Differential expression analysis with DESeq2 via tximport
// Calls scripts/deseq2_analysis.R

process DESEQ2 {
    label 'process_medium'

    conda     'envs/deseq2.yml'
    container 'quay.io/biocontainers/bioconductor-deseq2:1.42.0--r43hf17093f_2'

    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path quant_dirs     // collected Salmon quantification directories
    path samplesheet    // CSV with sample, condition, batch columns
    path tx2gene        // transcript-to-gene mapping (tx_id, gene_id, gene_name)

    output:
    path 'deseq2_results.tsv',       emit: results
    path 'deseq2_normalised_counts.tsv', emit: norm_counts
    path 'deseq2_size_factors.tsv',  emit: size_factors
    path 'plots/',                   emit: plots
    path 'deseq2_dds.rds',           emit: dds
    path '*.version.txt',            emit: version

    script:
    """
    Rscript ${projectDir}/scripts/deseq2_analysis.R \\
        --quant_dir    . \\
        --samplesheet  ${samplesheet} \\
        --tx2gene      ${tx2gene} \\
        --contrast     condition \\
        --numerator    dexamethasone \\
        --denominator  untreated \\
        --lfc_thresh   1 \\
        --padj_thresh  0.05 \\
        --threads      ${task.cpus} \\
        --outdir       .

    Rscript -e "cat(as.character(packageVersion('DESeq2')), '\\n')" > deseq2.version.txt
    """

    stub:
    """
    touch deseq2_results.tsv deseq2_normalised_counts.tsv deseq2_size_factors.tsv deseq2_dds.rds deseq2.version.txt
    mkdir -p plots
    """
}
