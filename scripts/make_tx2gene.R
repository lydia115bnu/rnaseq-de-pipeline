#!/usr/bin/env Rscript
# =============================================================
# make_tx2gene.R
# Build tx2gene.tsv from Ensembl GTF or biomaRt
# Run once before the pipeline:
#   Rscript scripts/make_tx2gene.R --gtf Homo_sapiens.GRCh38.110.gtf.gz --outfile data/tx2gene.tsv
# =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option("--gtf",     type="character", help="Path to Ensembl GTF file"),
  make_option("--outfile", type="character", default="data/tx2gene.tsv")
)
opt <- parse_args(OptionParser(option_list=option_list))

message("Parsing GTF: ", opt$gtf)
gtf <- rtracklayer::import(opt$gtf)
tx  <- as.data.frame(gtf) |>
  filter(type == "transcript") |>
  select(tx_id = transcript_id,
         gene_id,
         gene_name) |>
  distinct()

write_tsv(tx, opt$outfile, col_names=FALSE)
message("tx2gene written: ", opt$outfile, " (", nrow(tx), " transcripts)")
