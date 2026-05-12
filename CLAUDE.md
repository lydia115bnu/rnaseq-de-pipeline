# CLAUDE.md — rnaseq-de-pipeline

## What this pipeline does

End-to-end Nextflow DSL2 RNA-seq differential expression pipeline for the GSE52778 dataset (dexamethasone treatment of airway smooth muscle cells, 8 paired-end samples across 4 cell lines).

**Workflow:** FastQC → fastp → Salmon (index + quant) → MultiQC → DESeq2 → clusterProfiler (optional)

## Repository layout

```
rnaseq-de-pipeline/
├── main.nf                        # Nextflow entry point; full workflow definition
├── nextflow.config                # Profiles: docker | singularity | conda | test | slurm
├── conf/
│   ├── base.config                # Default resource labels (process_low/medium/high)
│   ├── slurm.config               # HPC SLURM executor settings
│   └── test.config                # Smoke-test overrides (chr22, 2 samples, no enrichment)
├── modules/
│   ├── fastqc.nf                  # Raw QC
│   ├── fastp.nf                   # Adapter trimming
│   ├── salmon.nf                  # SALMON_INDEX + SALMON_QUANT processes
│   ├── multiqc.nf                 # Aggregate QC report
│   ├── deseq2.nf                  # Calls scripts/deseq2_analysis.R
│   └── clusterprofiler.nf         # Calls scripts/pathway_enrichment.R
├── scripts/
│   ├── deseq2_analysis.R          # tximport → DESeq2 → PCA/volcano/heatmap/TSVs
│   ├── pathway_enrichment.R       # GO ORA, KEGG ORA, GSEA via clusterProfiler/fgsea
│   └── make_tx2gene.R             # Helper to build tx2gene.tsv from a GTF/Ensembl
├── envs/
│   ├── rnaseq.yml                 # conda: Salmon 1.10.3, FastQC, fastp, MultiQC
│   ├── deseq2.yml                 # conda: R 4.3.2, DESeq2, tximport, ggplot2, etc.
│   └── enrichment.yml             # conda: clusterProfiler, fgsea, enrichplot
├── data/
│   ├── samplesheet.csv            # 8-sample full run (SRR1039508–SRR1039517)
│   └── ids.csv                    # SRA accession list for download
└── results/                       # Output directories (gitkeep placeholders)
```

## Running the pipeline

### Prerequisites
- Nextflow ≥ 23.04.0
- One of: Docker, Singularity, or conda/mamba
- A human transcriptome FASTA (e.g., Ensembl `Homo_sapiens.GRCh38.cdna.all.fa.gz`) and genome FASTA (for decoy-aware Salmon index), OR a pre-built Salmon index

### Full run (conda profile example)
```bash
cd ~/Claude/rnaseq-de-pipeline

nextflow run main.nf \
  -profile conda \
  --transcriptome /path/to/Homo_sapiens.GRCh38.cdna.all.fa.gz \
  --genome        /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --input         data/samplesheet.csv \
  --tx2gene       data/tx2gene.tsv \
  --outdir        results
```

### Smoke test (no real data needed — requires test index & samplesheet_test.csv)
```bash
nextflow run main.nf -profile test,conda
```

### Skip enrichment
```bash
nextflow run main.nf -profile conda --skip_enrichment true ...
```

### Re-use a pre-built Salmon index
```bash
nextflow run main.nf -profile conda --salmon_index /path/to/index ...
```

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | `data/samplesheet.csv` | CSV with columns: sample, fastq_1, fastq_2, condition, batch |
| `--transcriptome` | null | cdna FASTA (required unless `--salmon_index` given) |
| `--genome` | null | Genome FASTA for decoy index (pass with `--transcriptome`) |
| `--salmon_index` | null | Pre-built index path; skips SALMON_INDEX step |
| `--tx2gene` | `data/tx2gene.tsv` | TSV: tx_id, gene_id, gene_name |
| `--organism` | `hsa` | KEGG organism code |
| `--org_db` | `org.Hs.eg.db` | Bioconductor OrgDb package |
| `--lfc_thresh` | `1.0` | |log2FC| cutoff for significance |
| `--padj_thresh` | `0.05` | FDR cutoff |
| `--skip_enrichment` | `false` | Set `true` to skip clusterProfiler step |

## Samplesheet format

```csv
sample,fastq_1,fastq_2,condition,batch
SRR1039508,data/reads/SRR1039508_1.fastq.gz,data/reads/SRR1039508_2.fastq.gz,untreated,N61311
SRR1039509,data/reads/SRR1039509_1.fastq.gz,data/reads/SRR1039509_2.fastq.gz,dexamethasone,N61311
```
- `condition` must match `--numerator` / `--denominator` in the DESeq2 call (defaults: dexamethasone vs untreated)
- `batch` is used as a covariate in the DESeq2 design (`~ batch + condition`)

## DESeq2 outputs (results/deseq2/)

| File | Description |
|---|---|
| `deseq2_results.tsv` | All genes: gene_id, gene_name, LFC (ashr-shrunk), padj, etc. |
| `deseq2_normalised_counts.tsv` | DESeq2 size-factor-normalised counts |
| `deseq2_size_factors.tsv` | Per-sample size factors |
| `deseq2_dds.rds` | Full DESeqDataSet R object |
| `plots/pca.pdf` | PCA coloured by condition + batch |
| `plots/volcano.pdf` | Volcano with top 20 DEG labels |
| `plots/heatmap_top50.pdf` | Heatmap of top 50 DEGs (VST) |
| `plots/pvalue_histogram.png` | Diagnostic p-value distribution |

## Enrichment outputs (results/enrichment/)

| File | Description |
|---|---|
| `go_enrichment.tsv` | GO BP ORA results (simplified) |
| `kegg_enrichment.tsv` | KEGG ORA results |
| `gsea_results.tsv` | GSEA GO BP results |
| `plots/go_ora_dotplot.pdf` | GO dotplot (top 20 terms) |
| `plots/go_ora_emap.pdf` | GO term similarity network |
| `plots/kegg_ora_dotplot.pdf` | KEGG dotplot |
| `plots/gsea_ridgeplot.pdf` | GSEA ridge plot |
| `plots/gsea_top_term.pdf` | GSEA enrichment plot for top term |

## Conda environments

Install all environments before a non-Docker/Singularity run:
```bash
mamba env create -f envs/rnaseq.yml
mamba env create -f envs/deseq2.yml
mamba env create -f envs/enrichment.yml
```

## Git / GitHub Desktop workflow

The repo lives at `~/Claude/rnaseq-de-pipeline`. After making changes:
1. Stage and commit via `git` CLI or GitHub Desktop.
2. Push to the remote (`RNAseq-pipeline` repo) so GitHub Desktop can pull.
3. Never force-push to `main`/`master`.

## Validation checklist

Before calling the pipeline production-ready, confirm:
- [ ] `nextflow run main.nf -profile test,conda` exits 0
- [ ] All 7 Nextflow modules parse without syntax errors (`nextflow inspect main.nf`)
- [ ] `data/tx2gene.tsv` exists (generate with `Rscript scripts/make_tx2gene.R`)
- [ ] `data/samplesheet_test.csv` and `data/test_index/` exist for the test profile
- [ ] `results/pipeline_info/execution_report.html` is generated after a run

## Missing data (not tracked in git)

- `data/reads/*.fastq.gz` — download from SRA using accessions in `data/ids.csv`
- `data/tx2gene.tsv` — generate with `scripts/make_tx2gene.R` from an Ensembl GTF
- `data/test_index/` — pre-built Salmon index for chr22 (for smoke tests)
- `data/samplesheet_test.csv` — 2-sample subset pointing to test reads
