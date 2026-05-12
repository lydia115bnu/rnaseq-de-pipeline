# rnaseq-de-pipeline

Nextflow DSL2 RNA-seq differential expression pipeline for the GSE52778 dataset (dexamethasone treatment of airway smooth muscle cells, 8 paired-end samples across 4 cell lines).

**Workflow:** FastQC → fastp → Salmon (index + quant) → MultiQC → DESeq2 → clusterProfiler (optional)

## Prerequisites

- Nextflow ≥ 23.04.0
- One of: Docker, Singularity, or conda/mamba

## Samplesheet format

`data/samplesheet.csv` — five columns, paired-end FASTQ:

```csv
sample,fastq_1,fastq_2,condition,batch
SRR1039508,data/reads/SRR1039508_1.fastq.gz,data/reads/SRR1039508_2.fastq.gz,untreated,N61311
SRR1039509,data/reads/SRR1039509_1.fastq.gz,data/reads/SRR1039509_2.fastq.gz,dexamethasone,N61311
```

`batch` is used as a covariate in the DESeq2 design (`~ batch + condition`).

## Usage

```bash
# Full run — build index from transcriptome FASTA
nextflow run main.nf \
  -profile conda \
  --transcriptome /path/to/Homo_sapiens.GRCh38.cdna.all.fa.gz \
  --genome        /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  --input         data/samplesheet.csv \
  --tx2gene       data/tx2gene.tsv \
  --outdir        results

# Reuse a pre-built Salmon index
nextflow run main.nf -profile conda --salmon_index /path/to/index ...

# Skip enrichment
nextflow run main.nf -profile conda --skip_enrichment true ...

# Smoke test (requires data/test_index and data/samplesheet_test.csv)
nextflow run main.nf -profile test,conda
```

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | `data/samplesheet.csv` | Samplesheet CSV |
| `--transcriptome` | required* | cdna FASTA for Salmon index |
| `--genome` | required* | Genome FASTA for decoy index |
| `--salmon_index` | null | Pre-built index (skips indexing step) |
| `--tx2gene` | `data/tx2gene.tsv` | TSV: tx_id, gene_id, gene_name |
| `--organism` | `hsa` | KEGG organism code |
| `--org_db` | `org.Hs.eg.db` | Bioconductor OrgDb package |
| `--lfc_thresh` | `1.0` | \|log2FC\| cutoff |
| `--padj_thresh` | `0.05` | FDR cutoff |
| `--skip_enrichment` | `false` | Skip clusterProfiler step |

\* Required unless `--salmon_index` is provided.

## Outputs

```
results/
  fastqc/                 FastQC reports (raw reads)
  fastp/                  Trimmed reads + QC reports
  salmon/                 Per-sample Salmon quant directories
  multiqc/                Aggregated MultiQC HTML report
  deseq2/
    deseq2_results.tsv              All genes: LFC (ashr-shrunk), padj
    deseq2_normalised_counts.tsv    Size-factor-normalised counts
    deseq2_size_factors.tsv         Per-sample size factors
    deseq2_dds.rds                  Full DESeqDataSet R object
    plots/pca.pdf                   PCA by condition + batch
    plots/volcano.pdf               Volcano with top DEG labels
    plots/heatmap_top50.pdf         Top 50 DEG heatmap (VST)
    plots/pvalue_histogram.png      Diagnostic p-value distribution
  enrichment/
    go_enrichment.tsv               GO BP ORA (simplified)
    kegg_enrichment.tsv             KEGG ORA
    gsea_results.tsv                GSEA GO BP
    plots/                          Dot plots, emap, ridge plots
  pipeline_info/                    Nextflow execution reports + DAG
```

## Conda environments

```bash
mamba env create -f envs/rnaseq.yml      # FastQC, fastp, Salmon, MultiQC
mamba env create -f envs/deseq2.yml      # DESeq2, tximport, ggplot2
mamba env create -f envs/enrichment.yml  # clusterProfiler, fgsea, org.Hs.eg.db
```

## Building tx2gene.tsv

Run once before the first pipeline execution:

```bash
Rscript scripts/make_tx2gene.R \
  --gtf     /path/to/Homo_sapiens.GRCh38.110.gtf.gz \
  --outfile data/tx2gene.tsv
```

## Reference genome downloads (Ensembl 110, GRCh38)

```bash
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
```
