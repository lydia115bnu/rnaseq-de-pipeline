// modules/clusterprofiler.nf
// Gene Ontology and KEGG pathway enrichment via clusterProfiler + fgsea

process ENRICHMENT {
    label 'process_medium'

    conda     'envs/enrichment.yml'
    container 'quay.io/biocontainers/bioconductor-clusterprofiler:4.10.0--r43hdfd78af_0'

    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path de_results      // deseq2_results.tsv with log2FC, padj, gene columns
    val  organism        // e.g. 'hsa' (KEGG) / 'Homo sapiens' (GO)
    val  org_db          // e.g. 'org.Hs.eg.db'

    output:
    path 'go_enrichment.tsv',        emit: go
    path 'kegg_enrichment.tsv',      emit: kegg
    path 'gsea_results.tsv',         emit: gsea
    path 'plots/',                   emit: plots
    path '*.version.txt',            emit: version

    script:
    """
    Rscript ${projectDir}/scripts/pathway_enrichment.R \\
        --de_results  ${de_results} \\
        --organism    ${organism} \\
        --org_db      ${org_db} \\
        --padj_thresh 0.05 \\
        --lfc_thresh  1 \\
        --outdir      .

    Rscript -e "cat(as.character(packageVersion('clusterProfiler')), '\\n')" > clusterprofiler.version.txt
    """

    stub:
    """
    touch go_enrichment.tsv kegg_enrichment.tsv gsea_results.tsv clusterprofiler.version.txt
    mkdir -p plots
    """
}
