#!/usr/bin/env Rscript
# =============================================================
# deseq2_analysis.R
# Differential expression analysis via tximport + DESeq2
# Called from modules/deseq2.nf
#
# Usage:
#   Rscript deseq2_analysis.R \
#     --quant_dir    <salmon_quant_dirs_parent> \
#     --samplesheet  <samplesheet.csv> \
#     --tx2gene      <tx2gene.tsv> \
#     --contrast     condition \
#     --numerator    dexamethasone \
#     --denominator  untreated \
#     --lfc_thresh   1 \
#     --padj_thresh  0.05 \
#     --threads      4 \
#     --outdir       .
# =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

# ── Argument parsing ──────────────────────────────────────────
option_list <- list(
  make_option("--quant_dir",   type="character", help="Directory containing per-sample Salmon output folders"),
  make_option("--samplesheet", type="character", help="CSV: sample, condition, batch columns"),
  make_option("--tx2gene",     type="character", help="TSV: tx_id, gene_id, gene_name"),
  make_option("--contrast",    type="character", default="condition"),
  make_option("--numerator",   type="character", help="Test condition (e.g. dexamethasone)"),
  make_option("--denominator", type="character", help="Reference condition (e.g. untreated)"),
  make_option("--lfc_thresh",  type="double",    default=1.0,  help="LFC threshold for results()"),
  make_option("--padj_thresh", type="double",    default=0.05, help="Adjusted p-value cutoff"),
  make_option("--threads",     type="integer",   default=4,    help="BiocParallel threads"),
  make_option("--outdir",      type="character", default=".",  help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(file.path(opt$outdir, "plots"), showWarnings=FALSE, recursive=TRUE)

# ── 1. Load samplesheet ──────────────────────────────────────
message("Loading samplesheet: ", opt$samplesheet)
meta <- read_csv(opt$samplesheet, show_col_types=FALSE)
meta$condition <- factor(meta$condition, levels=c(opt$denominator, opt$numerator))
meta$batch     <- factor(meta$batch)

# ── 2. Locate Salmon quant.sf files ─────────────────────────
quant_files <- file.path(opt$quant_dir, meta$sample, "quant.sf")
names(quant_files) <- meta$sample

missing <- quant_files[!file.exists(quant_files)]
if (length(missing) > 0) {
  stop("Missing quant.sf files:\n", paste(missing, collapse="\n"))
}

# ── 3. tx2gene mapping ───────────────────────────────────────
message("Loading tx2gene: ", opt$tx2gene)
tx2gene <- read_tsv(opt$tx2gene, col_names=c("tx_id","gene_id","gene_name"),
                    show_col_types=FALSE)

# ── 4. tximport: import Salmon output ────────────────────────
message("Importing Salmon quantifications with tximport...")
txi <- tximport(
  quant_files,
  type                 = "salmon",
  tx2gene              = tx2gene[, c("tx_id","gene_id")],
  countsFromAbundance  = "lengthScaledTPM",
  ignoreAfterBar       = TRUE  # handles Ensembl versioned IDs
)

# ── 5. Build DESeqDataSet ────────────────────────────────────
message("Building DESeqDataSet...")
dds <- DESeqDataSetFromTximport(
  txi,
  colData = meta,
  design  = ~ batch + condition   # account for cell-line batch
)

# Pre-filter: keep genes with >=10 counts in at least 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message(sprintf("Genes retained after filtering: %d / %d", sum(keep), length(keep)))

# ── 6. Run DESeq2 ────────────────────────────────────────────
message("Running DESeq2...")
BiocParallel::register(BiocParallel::MulticoreParam(opt$threads))
dds <- DESeq(dds, parallel=TRUE)

# Check size factors — warn if extreme
sf <- sizeFactors(dds)
message("Size factors range: ", round(min(sf),2), " – ", round(max(sf),2))
if (max(sf)/min(sf) > 5) {
  warning("Size factor range >5× — check for library preparation outliers")
}

# ── 7. Extract results with LFC shrinkage ────────────────────
message("Extracting results: ", opt$numerator, " vs ", opt$denominator)
res_raw <- results(
  dds,
  contrast     = c(opt$contrast, opt$numerator, opt$denominator),
  lfcThreshold = 0,         # use raw LFC for shrinkage input
  alpha        = opt$padj_thresh
)

# LFC shrinkage with ashr (recommended over apeglm for contrast-based)
res_shrunk <- lfcShrink(
  dds,
  contrast = c(opt$contrast, opt$numerator, opt$denominator),
  res      = res_raw,
  type     = "ashr"
)

# Add gene symbols from rownames (assuming ENSEMBL IDs as rownames)
res_df <- as.data.frame(res_shrunk) |>
  tibble::rownames_to_column("gene_id") |>
  left_join(tx2gene |> distinct(gene_id, gene_name), by="gene_id") |>
  arrange(padj)

# Diagnostic: p-value histogram
png(file.path(opt$outdir, "plots", "pvalue_histogram.png"), width=800, height=600)
hist(res_raw$pvalue, breaks=50, col="steelblue", border="white",
     main=paste0("p-value distribution: ", opt$numerator, " vs ", opt$denominator),
     xlab="Raw p-value")
dev.off()

# ── 8. Significant DEGs ──────────────────────────────────────
sig <- res_df |>
  filter(!is.na(padj),
         padj < opt$padj_thresh,
         abs(log2FoldChange) >= opt$lfc_thresh)

message(sprintf("Significant DEGs (padj<%.2f, |LFC|>=%.1f): %d",
                opt$padj_thresh, opt$lfc_thresh, nrow(sig)))
message(sprintf("  Upregulated  : %d", sum(sig$log2FoldChange > 0)))
message(sprintf("  Downregulated: %d", sum(sig$log2FoldChange < 0)))

# ── 9. PCA plot ──────────────────────────────────────────────
vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup=c("condition","batch"), returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color=condition, shape=batch, label=name)) +
  geom_point(size=4) +
  geom_text_repel(size=3, show.legend=FALSE) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values=c("untreated"="#4393c3","dexamethasone"="#d6604d")) +
  theme_bw(14) +
  ggtitle("PCA — VST-normalised counts")

ggsave(file.path(opt$outdir, "plots", "pca.pdf"), p_pca, width=7, height=6)

# ── 10. Volcano plot ─────────────────────────────────────────
volcano_df <- res_df |>
  mutate(
    significance = case_when(
      padj < opt$padj_thresh & log2FoldChange >=  opt$lfc_thresh ~ "Up",
      padj < opt$padj_thresh & log2FoldChange <= -opt$lfc_thresh ~ "Down",
      TRUE ~ "NS"
    ),
    label = ifelse(significance != "NS" &
                   rank(padj, na.last="keep") <= 20, gene_name, NA)
  )

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, -log10(padj), color=significance)) +
  geom_point(alpha=0.6, size=1.2) +
  geom_text_repel(aes(label=label), size=2.8, max.overlaps=15) +
  scale_color_manual(values=c("Up"="#d6604d","Down"="#4393c3","NS"="grey70")) +
  geom_hline(yintercept=-log10(opt$padj_thresh), linetype="dashed", color="grey40") +
  geom_vline(xintercept=c(-opt$lfc_thresh, opt$lfc_thresh), linetype="dashed", color="grey40") +
  theme_bw(14) +
  labs(title=paste0(opt$numerator, " vs ", opt$denominator),
       x="Shrunken log2 fold change", y="-log10(adjusted p-value)", color="")

ggsave(file.path(opt$outdir, "plots", "volcano.pdf"), p_volcano, width=8, height=7)

# ── 11. Heatmap of top DEGs ──────────────────────────────────
top_genes <- head(sig$gene_id, 50)
if (length(top_genes) >= 2) {
  mat <- assay(vsd)[top_genes, ]
  rownames(mat) <- sig$gene_name[match(top_genes, sig$gene_id)]
  mat <- mat - rowMeans(mat)

  ann_col <- as.data.frame(colData(vsd)[, c("condition","batch")])

  pdf(file.path(opt$outdir, "plots", "heatmap_top50.pdf"), width=9, height=12)
  pheatmap(mat,
           annotation_col = ann_col,
           color          = colorRampPalette(rev(brewer.pal(9,"RdBu")))(100),
           show_rownames  = TRUE,
           show_colnames  = TRUE,
           fontsize_row   = 7,
           main           = paste0("Top 50 DEGs: ", opt$numerator, " vs ", opt$denominator))
  dev.off()
}

# ── 12. Save outputs ─────────────────────────────────────────
write_tsv(res_df,                                  file.path(opt$outdir, "deseq2_results.tsv"))
write_tsv(as.data.frame(counts(dds, normalized=TRUE)) |>
            tibble::rownames_to_column("gene_id"), file.path(opt$outdir, "deseq2_normalised_counts.tsv"))
write_tsv(tibble::tibble(sample=names(sizeFactors(dds)), size_factor=sizeFactors(dds)),
                                                   file.path(opt$outdir, "deseq2_size_factors.tsv"))
saveRDS(dds,                                       file.path(opt$outdir, "deseq2_dds.rds"))

message("DESeq2 analysis complete. Results written to: ", opt$outdir)
