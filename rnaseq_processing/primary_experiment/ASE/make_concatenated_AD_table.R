
suppressPackageStartupMessages({
  library(data.table)
})



allele_dir <- "/rs/rs_grp_gmccol/GMCcol_ASE_analysis/per_individual_joint_pileup/keep_zero_ALT_sites_analysis/ASE_input"
out_prefix <- "all_individuals"

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%F %T")))
  cat(sprintf(...), "\n")
}

files <- list.files(
  allele_dir,
  pattern = "_allele_counts_long\\.txt\\.gz$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No *_allele_counts_long.txt.gz files found in: ", allele_dir)
}

msg("Found %d allele depth files", length(files))

dt_list <- vector("list", length(files))

for (i in seq_along(files)) {
  if (i %% 25 == 0) msg("Reading %d / %d", i, length(files))
  dt_list[[i]] <- fread(files[i])
}

dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

if (!("Individual" %in% names(dt))) {
  dt[, Individual := sub("_allele_counts_long\\.txt\\.gz$", "", basename(files[1]))]
  warning("Individual column was missing. Check per-individual script.")
}

out_file <- file.path(
  allele_dir,
  paste0(out_prefix, "_allele_counts_long.txt.gz")
)

msg("Writing combined allele depth table: %s", out_file)

fwrite(
  dt,
  out_file,
  sep = "\t",
  compress = "gzip"
)

