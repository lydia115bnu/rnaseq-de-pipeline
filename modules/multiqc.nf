// modules/multiqc.nf
// Aggregate QC reports from FastQC, fastp, and Salmon

process MULTIQC {
    label 'process_low'

    conda     'envs/rnaseq.yml'
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files  // collected mix of FastQC zips, fastp JSONs, Salmon logs

    output:
    path 'multiqc_report.html',   emit: report
    path 'multiqc_data/',         emit: data
    path '*.version.txt',         emit: version

    script:
    """
    multiqc \\
        --force \\
        --title "GSE52778 Airway RNA-seq QC" \\
        --comment "Dexamethasone treatment of airway smooth muscle cells" \\
        .

    multiqc --version | sed 's/multiqc, version //' > multiqc.version.txt
    """

    stub:
    """
    touch multiqc_report.html multiqc.version.txt
    mkdir -p multiqc_data
    """
}
