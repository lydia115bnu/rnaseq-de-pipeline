#!/usr/bin/env nextflow
// ============================================================
// GSE52778 Airway RNA-seq Pipeline
// Dexamethasone treatment of airway smooth muscle cells
// ============================================================
nextflow.enable.dsl = 2

include { FASTQC        } from './modules/fastqc'
include { FASTP         } from './modules/fastp'
include { SALMON_INDEX  } from './modules/salmon'
include { SALMON_QUANT  } from './modules/salmon'
include { MULTIQC       } from './modules/multiqc'
include { DESEQ2        } from './modules/deseq2'
include { ENRICHMENT    } from './modules/clusterprofiler'

// ── Parameter defaults ─────────────────────────────────────
params.input          = 'data/samplesheet.csv'
params.outdir         = 'results'
params.transcriptome  = null          // path to cdna FASTA (required)
params.genome         = null          // path to genome FASTA (for decoy index)
params.salmon_index   = null          // pre-built index — skips SALMON_INDEX if provided
params.tx2gene        = 'data/tx2gene.tsv'
params.organism       = 'hsa'         // KEGG organism code
params.org_db         = 'org.Hs.eg.db'
params.skip_enrichment = false

// Thresholds
params.lfc_thresh     = 1.0
params.padj_thresh    = 0.05

// Validate required params
if (!params.transcriptome && !params.salmon_index) {
    error "Provide --transcriptome (FASTA) or --salmon_index (pre-built directory)"
}

// ── Helper: parse samplesheet ──────────────────────────────
def parse_samplesheet(csv_file) {
    Channel
        .fromPath(csv_file, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:        row.sample,
                condition: row.condition,
                batch:     row.batch
            ]
            def reads = [
                file(row.fastq_1, checkIfExists: true),
                file(row.fastq_2, checkIfExists: true)
            ]
            return [meta, reads]
        }
}

// ── Main workflow ──────────────────────────────────────────
workflow {

    // 1. Input channel
    reads_ch = parse_samplesheet(params.input)

    // 2. Raw QC
    FASTQC(reads_ch)

    // 3. Adapter trimming
    FASTP(reads_ch)

    // 4. Build Salmon index (or reuse existing)
    if (params.salmon_index) {
        salmon_index = file(params.salmon_index)
    } else {
        SALMON_INDEX(
            file(params.transcriptome, checkIfExists: true),
            file(params.genome,        checkIfExists: true)
        )
        salmon_index = SALMON_INDEX.out.index
    }

    // 5. Quantification — run on trimmed reads
    SALMON_QUANT(FASTP.out.trimmed, salmon_index)

    // 6. Aggregate QC
    qc_mix = FASTQC.out.zip
        .mix(FASTP.out.json)
        .mix(SALMON_QUANT.out.log)
        .map { meta, files -> files }
        .collect()

    MULTIQC(qc_mix)

    // 7. Differential expression
    DESEQ2(
        SALMON_QUANT.out.quants.map { meta, dir -> dir }.collect(),
        file(params.input,   checkIfExists: true),
        file(params.tx2gene, checkIfExists: true)
    )

    // 8. Pathway enrichment (optional)
    if (!params.skip_enrichment) {
        ENRICHMENT(
            DESEQ2.out.results,
            params.organism,
            params.org_db
        )
    }
}

// ── Completion handler ─────────────────────────────────────
workflow.onComplete {
    log.info """
    ╔══════════════════════════════════════╗
    ║   GSE52778 RNA-seq pipeline done!    ║
    ╠══════════════════════════════════════╣
    ║  Status  : ${workflow.success ? 'SUCCESS ✓' : 'FAILED  ✗'}
    ║  Duration: ${workflow.duration}
    ║  Results : ${params.outdir}
    ╚══════════════════════════════════════╝
    """.stripIndent()
}
