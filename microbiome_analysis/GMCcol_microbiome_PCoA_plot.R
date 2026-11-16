# PCoA of microbiome sample metagenomes using the covariate file with CLR taxa
# Color = Urban/Rural microbiome, Shape = Country
# No phyloseq needed: uses ape::pcoa() directly

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ape)
})

# --- paths ---
microbiome_dir <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/microbiome_analysis"
cov_dir <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates"

# --- Load covariates with CLR abundances for all 136 taxa ---
cv_micro <- read.table(
  file.path(cov_dir, "GMCcol_deep_covariates_all_microbial_taxa_7_29_2026.txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = ""
)

# Same non-taxon column list as your cASE script
non_taxa_cols <- c("Microbiome", "SampleName", "Urbanism",
                   "Treatment", "Colonoid_line", "Sample_ID", "SampleID", "Sorted",
                   "RawReads", "Index", "Quality", "Clean", "i7.index",
                   "i5.index", "Batch", "Colonoid_age", "Colonoid_Sex",
                   "Colonoid_Race", "Country", "Lifestyle", "Ratio", "TotalReads",
                   "NonDuplicated_reads", "Firmicutes", "Prevotella", "Bacteroides", "Shannon")

taxa_cols <- setdiff(colnames(cv_micro), non_taxa_cols)
cat("Number of taxa columns:", length(taxa_cols), "\n")   # expect 136

# --- Collapse to one row per microbiome sample ---
# (many rows per sample in the covariate file; taxa values are identical within a sample)
micro_samples <- cv_micro %>%
  filter(!is.na(SampleName)) %>%
  select(Microbiome, SampleName, Country, Urbanism, all_of(taxa_cols)) %>%
  distinct(SampleName, .keep_all = TRUE)

cat("Number of microbiome samples:", nrow(micro_samples), "\n")   # expect 8
stopifnot(nrow(micro_samples) == 8)

# Check that taxa values really are constant within each sample
check_const <- cv_micro %>%
  group_by(SampleName) %>%
  summarise(across(all_of(taxa_cols), ~ n_distinct(.x) == 1), .groups = "drop")
stopifnot(all(unlist(check_const[, taxa_cols])))

# --- Samples x taxa matrix ---
clr_mat <- micro_samples %>%
  select(SampleName, all_of(taxa_cols)) %>%
  column_to_rownames("SampleName") %>%
  as.matrix()
storage.mode(clr_mat) <- "numeric"

meta_pcoa <- micro_samples %>% select(Microbiome, SampleName, Country, Urbanism)

# --- Distance + PCoA ---
d <- dist(clr_mat, method = "euclidean")   # Aitchison distance on CLR data
pcoa_res <- pcoa(d)
rel_eig <- pcoa_res$values$Relative_eig * 100

# --- Sample coordinates + metadata ---
mycolors <- c("Urban Microbiome" = "#0072B2", "Rural Microbiome" = "#009E73")

pcoa_df <- pcoa_res$vectors[, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column("SampleName") %>%
  left_join(meta_pcoa, by = "SampleName") %>%
  mutate(Urbanism = factor(paste(Urbanism, "Microbiome"),
                           levels = c("Rural Microbiome", "Urban Microbiome")))

stopifnot(!any(is.na(pcoa_df$Country)), !any(is.na(pcoa_df$Urbanism)))

# --- Plot ---
n_country <- length(unique(pcoa_df$Country))

txt <- 16   # one size for all text

p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2,
                              fill = Urbanism, shape = Country)) +
  geom_point(size = 5, stroke = 0.75, color = "white") +
  scale_fill_manual(values = mycolors, name = NULL) +
  scale_shape_manual(values = c(21, 22, 24, 23, 25)[seq_len(n_country)],
                     name = NULL) +
  guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, size = 5)),
         shape = guide_legend(order = 2, override.aes = list(fill = "grey60", size = 5))) +
  xlab(sprintf("PC1 [%.1f%%]", rel_eig[1])) +
  ylab(sprintf("PC2 [%.1f%%]", rel_eig[2])) +
  ggtitle("PCoA of Microbiome Sample Metagenomes") +
  theme_classic(base_size = txt) +
  theme(
    text         = element_text(size = txt),
    plot.title   = element_text(size = txt, face = "bold", hjust = 0.5),
    plot.title.position = "plot",          # center over the whole plot, not just the panel
    axis.title   = element_text(size = txt),
    axis.text    = element_text(size = txt, color = "black"),
    legend.text  = element_text(size = txt),
    legend.position = "right"
  )

# With legend
ggsave(file.path(microbiome_dir, "PCoA_microbiome_samples_country_urbanism.pdf"),
       p_pcoa, width = 9, height = 5.5)

# Without legend (narrower so the panel stays about the same size)
pcoa2 <- p_pcoa + theme(legend.position = "none")
ggsave(file.path(microbiome_dir, "PCoA_microbiome_samples_country_urbanism_no_legend.pdf"),
       pcoa2, width = 6.5, height = 5.5)