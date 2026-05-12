#!/usr/bin/env Rscript
# =============================================================
# pathway_enrichment.R
# GO + KEGG ORA and GSEA using clusterProfiler & fgsea
# Called from modules/clusterprofiler.nf
#
# Usage:
#   Rscript pathway_enrichment.R \
#     --de_results  deseq2_results.tsv \
#     --organism    hsa \
#     --org_db      org.Hs.eg.db \
#     --padj_thresh 0.05 \
#     --lfc_thresh  1 \
#     --outdir      .
# =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# ── Argument parsing ──────────────────────────────────────────
option_list <- list(
  make_option("--de_results",  type="character", help="deseq2_results.tsv"),
  make_option("--organism",    type="character", default="hsa",          help="KEGG organism code"),
  make_option("--org_db",      type="character", default="org.Hs.eg.db", help="Bioc OrgDb package"),
  make_option("--padj_thresh", type="double",    default=0.05),
  make_option("--lfc_thresh",  type="double",    default=1.0),
  make_option("--outdir",      type="character", default=".")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir.create(file.path(opt$outdir, "plots"), showWarnings=FALSE, recursive=TRUE)

# Dynamically load the requested OrgDb
library(opt$org_db, character.only=TRUE)
org <- get(opt$org_db)

# ── 1. Load DE results ────────────────────────────────────────
message("Loading DE results: ", opt$de_results)
de <- read_tsv(opt$de_results, show_col_types=FALSE)

# Significant gene sets
sig_up   <- de |> filter(padj < opt$padj_thresh, log2FoldChange >=  opt$lfc_thresh)
sig_down <- de |> filter(padj < opt$padj_thresh, log2FoldChange <= -opt$lfc_thresh)
sig_all  <- bind_rows(sig_up, sig_down)

# Ranked gene list for GSEA (all genes, ranked by -log10(padj) * sign(LFC))
ranked <- de |>
  filter(!is.na(padj), !is.na(log2FoldChange)) |>
  mutate(rank_score = -log10(padj + 1e-300) * sign(log2FoldChange)) |>
  arrange(desc(rank_score))

# ── 2. Gene ID conversion (gene_name → ENTREZID) ─────────────
message("Converting gene IDs to ENTREZID...")
convert_ids <- function(genes, key="SYMBOL") {
  bitr(genes, fromType=key, toType="ENTREZID", OrgDb=org, drop=TRUE)
}

sig_entrez <- convert_ids(sig_all$gene_name)
ranked_entrez <- convert_ids(ranked$gene_name)

ranked_vec <- ranked |>
  inner_join(ranked_entrez, by=c("gene_name"="SYMBOL")) |>
  distinct(ENTREZID, .keep_all=TRUE) |>
  arrange(desc(rank_score))
gene_list <- setNames(ranked_vec$rank_score, ranked_vec$ENTREZID)

# ── 3. GO ORA (over-representation) ──────────────────────────
message("Running GO over-representation analysis...")
go_ora <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = ranked_entrez$ENTREZID,
  OrgDb         = org,
  ont           = "BP",     # Biological Process; also try "MF", "CC"
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

if (!is.null(go_ora) && nrow(go_ora) > 0) {
  # Simplify redundant GO terms
  go_ora_simplified <- simplify(go_ora, cutoff=0.7, by="p.adjust")

  # Dotplot
  p_go_dot <- dotplot(go_ora_simplified, showCategory=20) +
    ggtitle("GO Biological Process — ORA") +
    theme_bw(12)
  ggsave(file.path(opt$outdir, "plots", "go_ora_dotplot.pdf"), p_go_dot, width=10, height=8)

  # Emap
  go_ora_simplified <- pairwise_termsim(go_ora_simplified)
  p_emap <- emapplot(go_ora_simplified, showCategory=30) +
    ggtitle("GO BP — term similarity network")
  ggsave(file.path(opt$outdir, "plots", "go_ora_emap.pdf"), p_emap, width=12, height=10)

  write_tsv(as.data.frame(go_ora_simplified), file.path(opt$outdir, "go_enrichment.tsv"))
  message(sprintf("GO ORA: %d significant terms", nrow(go_ora_simplified)))
} else {
  message("No significant GO terms found")
  write_tsv(tibble::tibble(), file.path(opt$outdir, "go_enrichment.tsv"))
}

# ── 4. KEGG ORA ──────────────────────────────────────────────
message("Running KEGG over-representation analysis...")
kegg_ora <- enrichKEGG(
  gene          = sig_entrez$ENTREZID,
  universe      = ranked_entrez$ENTREZID,
  organism      = opt$organism,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

if (!is.null(kegg_ora) && nrow(kegg_ora) > 0) {
  p_kegg <- dotplot(kegg_ora, showCategory=20) +
    ggtitle("KEGG pathway — ORA") +
    theme_bw(12)
  ggsave(file.path(opt$outdir, "plots", "kegg_ora_dotplot.pdf"), p_kegg, width=10, height=8)

  write_tsv(as.data.frame(kegg_ora), file.path(opt$outdir, "kegg_enrichment.tsv"))
  message(sprintf("KEGG ORA: %d significant pathways", nrow(kegg_ora)))
} else {
  message("No significant KEGG pathways found")
  write_tsv(tibble::tibble(), file.path(opt$outdir, "kegg_enrichment.tsv"))
}

# ── 5. GSEA (Gene Set Enrichment Analysis) ───────────────────
message("Running GSEA...")
gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org,
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE,
  eps           = 0
)

if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  gsea_go <- setReadable(gsea_go, OrgDb=org, keyType="ENTREZID")

  # Ridge plot of top enriched sets
  p_ridge <- ridgeplot(gsea_go, showCategory=20) +
    ggtitle("GSEA GO BP — enrichment distribution") +
    theme_bw(12)
  ggsave(file.path(opt$outdir, "plots", "gsea_ridgeplot.pdf"), p_ridge, width=12, height=10)

  # GSEA plot for top term
  top_term <- gsea_go$ID[1]
  p_gsea1  <- gseaplot2(gsea_go, geneSetID=top_term, title=gsea_go$Description[1])
  ggsave(file.path(opt$outdir, "plots", "gsea_top_term.pdf"), p_gsea1, width=9, height=7)

  write_tsv(as.data.frame(gsea_go), file.path(opt$outdir, "gsea_results.tsv"))
  message(sprintf("GSEA: %d significant gene sets", nrow(gsea_go)))
} else {
  message("No significant GSEA gene sets found")
  write_tsv(tibble::tibble(), file.path(opt$outdir, "gsea_results.tsv"))
}

message("Pathway enrichment complete. Results written to: ", opt$outdir)
