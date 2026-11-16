#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

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

num_samples <- length(samples)
num_SNPs <- nrow(dt)

# Pivot longer
long <- melt(
  dt,
  id.vars = c("chr", "pos", "rsID", "ref", "alt"),
  measure.vars = samples,
  variable.name = "Sample_ID",
  value.name = "AD"
)

long[, ref_count := parse_ref(AD)]
long[, alt_count := parse_alt(AD)]
long[, total := ref_count + alt_count]
long[, alt_ratio := fifelse(total > 0, alt_count / total, NA_real_)]

long[, Individual := ind]

setcolorder(
  long,
  c("Individual", "chr", "pos", "rsID", "ref", "alt",
    "Sample_ID", "AD", "ref_count", "alt_count", "total", "alt_ratio")
)

# QC summary across all SNP:sample pairs
qc_summary <- data.table(
  Individual = ind,
  num_samples = num_samples,
  num_SNPs = num_SNPs,
  num_pairs = nrow(long),
  ref_only_pairs = sum(long$ref_count > 0 & long$alt_count == 0),
  alt_only_pairs = sum(long$ref_count == 0 & long$alt_count > 0),
  both_pairs = sum(long$ref_count > 0 & long$alt_count > 0),
  no_reads_pairs = sum(long$ref_count == 0 & long$alt_count == 0),
  pairs_total_gt0 = sum(long$total > 0)
)

# Save outputs
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fwrite(
  long,
  file.path(out_dir, paste0(ind, "_allele_counts_long.txt.gz")),
  sep = "\t",
  compress = "gzip"
)

fwrite(
  qc_summary,
  file.path(out_dir, paste0(ind, "_QC_summary.tsv")),
  sep = "\t"
)

msg("Finished allele depth table for %s", ind)
msg("SNPs: %d | Samples: %d | SNP:sample pairs: %d",
    num_SNPs, num_samples, nrow(long))
msg("Pairs with reads: %d | no reads: %d",
    qc_summary$pairs_total_gt0, qc_summary$no_reads_pairs)