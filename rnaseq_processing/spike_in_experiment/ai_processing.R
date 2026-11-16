#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: ai_processing.R <IND> <pileup_dir> <out_dir>")
}

ind <- args[1]
pileup_dir <- args[2]
out_dir <- args[3]

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%F %T")))
  cat(sprintf(...), "\n")
}


# Load data
ad_wide <- file.path(pileup_dir, paste0(ind, ".AD.wide.tsv.gz"))
samples <- fread(file.path(pileup_dir, paste0(ind, ".samples.txt")), header = FALSE)$V1

dt <- fread(
  ad_wide,
  header = FALSE,
  colClasses = rep("character", 5 + length(samples))
)

setnames(dt, c("chr", "pos", "rsID", "ref", "alt", samples))


# AD parsing
parse_ref <- function(x) {
  x[is.na(x)] <- "0"
  out <- as.integer(sub(",.*$", "", x))
  out[is.na(out)] <- 0L
  out
}

parse_alt <- function(x) {
  x[is.na(x)] <- "0"
  has_comma <- grepl(",", x, fixed = TRUE)
  out <- integer(length(x))
  out[has_comma] <- as.integer(sub("^.*?,", "", x[has_comma]))
  out[is.na(out)] <- 0L
  out
}

sample_cols <- samples
num_samples <- length(sample_cols)
num_SNPs_before_filter <- nrow(dt)

ad_mat <- as.matrix(dt[, ..sample_cols])

ref_counts <- matrix(
  parse_ref(as.vector(ad_mat)),
  nrow = nrow(dt),
  ncol = num_samples
)

alt_counts <- matrix(
  parse_alt(as.vector(ad_mat)),
  nrow = nrow(dt),
  ncol = num_samples
)


# Pivot longer
long <- melt(
  dt,
  id.vars = c("chr", "pos", "rsID", "ref", "alt"),
  measure.vars = sample_cols,
  variable.name = "Sample_ID",
  value.name = "AD"
)

long[, ref_count := parse_ref(AD)]
long[, alt_count := parse_alt(AD)]
long[, total := ref_count + alt_count]
long[, alt_ratio := fifelse(total > 0, alt_count / total, NA_real_)]


# QC summary
qc_summary <- data.table(
  Individual = ind,
  num_samples = num_samples,
  num_SNPs_before_filter = num_SNPs_before_filter,
  ref_only_pairs = sum(ref_counts > 0 & alt_counts == 0),
  alt_only_pairs = sum(ref_counts == 0 & alt_counts > 0),
  both_pairs = sum(ref_counts > 0 & alt_counts > 0),
  no_reads_pairs = sum(ref_counts == 0 & alt_counts == 0)
)


# Save outputs
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(
  long,
  file.path(out_dir, paste0(ind, "_allele_counts_long.txt.gz")),
  compress = "gzip"
)

fwrite(
  qc_summary,
  file.path(out_dir, paste0(ind, "_QC_summary.tsv")),
  sep = "\t"
)

msg("Finished ASE processing for %s", ind)
msg("SNPs before filter: %d", num_SNPs_before_filter)
msg("Rows: %d | Samples represented: %d", nrow(long), uniqueN(long$Sample_ID))