# Script to generate covariate file with all microbial abundances and heatmap
# GMCcol Project


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)



# --- paths ---
microbiome_dir <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/microbiome_analysis"
cov_dir <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates"

# --- Load data ---
clr <- read.csv(
    "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/microbiome_analysis/clr_species.csv",
    row.names=1,
    check.names=FALSE
)

tax <- read_tsv("/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/microbiome_analysis/taxonomic_names_updates.tsv")

tax <- tax %>%
    mutate(oldTaxonomy=gsub("\\.", "", oldTaxonomy),
        new_name=case_when(
        species != "" & species != "s__" ~ gsub("s__", "", species),
        TRUE ~ gsub("g__", "", genus)))

name_map <- setNames(tax$new_name, tax$oldTaxonomy)

rownames(clr) <- sapply(rownames(clr), function(taxon) {
    if (taxon %in% names(name_map)) name_map[[taxon]] else taxon
})


# --- Load covariates ---
cv <- read.table(
  file.path(cov_dir, "GMCcol_deep_covariates.txt"),
  sep = "\t",
  stringsAsFactors = FALSE,
  header = TRUE,
  comment = ""
)

# --- Sample mapping ---
SampleName <- c(
  "control", "6199II", "4041OU", "5905ZO", "5744TZ",
  "9261PK", "8283LS", "9794IK", "6310FW"
)

Microbiome <- c(
  "control", "GM29", "GM11", "GM43", "GM41",
  "GM06", "GM04", "GM45", "GM33"
)

microbiome_sample_ids <- data.frame(SampleName, Microbiome)

cv <- merge(cv, microbiome_sample_ids, by = "Microbiome")

cv_filt <- cv %>%
  filter(Treatment == "microbiome") %>%
  select(Microbiome, SampleName, Urbanism) %>%
  unique()

# --- Subset CLR table to microbiome samples only ---
my_samples <- cv_filt$SampleName
clr_sub <- clr[, my_samples]

cat("Number of taxa:", nrow(clr_sub), "\n")
# Number of taxa: 136
cat("Number of microbiome samples:", ncol(clr_sub), "\n")
# Number of microbiome samples: 8


# --- Spearman correlation matrix across taxa ---
cor_mat_all <- cor(t(clr_sub), method = "spearman")

# --- Color palette ---
cor_cols <- colorRampPalette(c("blue", "white", "red"))(100)
cor_breaks <- seq(-1, 1, length.out = 101)


# ============================================================
# Explicit clustering by Spearman correlation
# Distance = 1 - rho
# Check which linkage algorithm has highest correlation
# ============================================================

dist_signed <- as.dist(1 - cor_mat_all)
hc_signed_ward <- hclust(dist_signed, method = "ward.D2")
hc_signed_avg <- hclust(dist_signed, method = "average")
hc_complete <- hclust(dist_signed, method = "complete")

cor(dist_signed, cophenetic(hc_signed_ward))
#[1] 0.7745268
cor(dist_signed, cophenetic(hc_signed_avg))
#[1] 0.8015153
cor(dist_signed, cophenetic(hc_complete))
#[1] 0.67352


# ============================================================
# Explicit clustering by Spearman correlation
# Distance = 1 - rho
# Positive correlations cluster together; negative correlations separate.
# use average linkage as it has highest copenhetic correlation
# ============================================================

pdf(
  file.path(microbiome_dir, "taxa_spearman_correlation_heatmap_dist_1minusR_average.pdf"),
  width = 22,
  height = 20
)

p_signed <- pheatmap(
  cor_mat_all,
  color = cor_cols,
  breaks = cor_breaks,
  cluster_rows = hc_signed_avg,
  cluster_cols = hc_signed_avg,
  fontsize_row = 4,
  fontsize_col = 4,
  main = "Taxa Spearman correlation heatmap | distance = 1 - rho | Average"
)

print(p_signed)

dev.off()


# ============================================================
# Candidate cluster cuts using distance = 1 - rho
# Use these to visually decide the cluster number.
# ============================================================

#candidate_k <- c(5, 10, 15, 20, 25, 30, 35, 40)

candidate_k <- c(2,3, 4, 6, 7, 8, 9, 10, 11 , 12, 13, 14, 15)

for (k in candidate_k) {
  
  pdf(
    file.path(
      microbiome_dir,
      paste0("taxa_spearman_correlation_heatmap_dist_1minusR_average_k", k, ".pdf")
    ),
    width = 22,
    height = 20
  )
  
  pheatmap(
    cor_mat_all,
    color = cor_cols,
    breaks = cor_breaks,
    cluster_rows = hc_signed_avg,
    cluster_cols = hc_signed_avg,
    cutree_rows = k,
    cutree_cols = k,
    fontsize_row = 4,
    fontsize_col = 4,
    main = paste0(
      "Taxa Spearman correlation heatmap | distance = 1 - rho | Average Linkage | k = ",
      k
    )
  )
  
  dev.off()
}


# ============================================================
# Save dendrogram objects
# ============================================================

saveRDS(
  hc_signed_avg,
  file.path(microbiome_dir, "taxa_spearman_hclust_dist_1minusR_average.rds")
)


library(cluster)
sil_widths <- sapply(2:20, function(k) {
  clusters <- cutree(hc_signed_avg, k = k)
  mean(silhouette(clusters, dist_signed)[, 3])
})


pdf(
  file.path(microbiome_dir, "silhouette_width_by_k_avg.pdf"),
  width = 7,
  height = 5
)
plot(2:20, sil_widths, type = "b", xlab = "k", ylab = "mean silhouette width",
     main = "Silhouette width vs. number of clusters")
dev.off()
#peaks again at 12

cluster_labels_12 <- cutree(hc_signed_avg, k = 12)

tightness_12 <- sapply(sort(unique(cluster_labels_12)), function(cl) {
  taxa_in_cl <- names(cluster_labels_12)[cluster_labels_12 == cl]
  if (length(taxa_in_cl) < 2) return(NA)
  submat <- cor_mat_all[taxa_in_cl, taxa_in_cl]
  mean(submat[upper.tri(submat)])
})

n_taxa_12 <- table(cluster_labels_12)

tightness_df <- data.frame(
  cluster = sort(unique(cluster_labels_12)),
  n_taxa = as.integer(n_taxa_12[as.character(sort(unique(cluster_labels_12)))]),
  mean_within_rho = tightness_12
)




cluster_df_12 <- data.frame(
  taxon = names(cluster_labels_12),
  cluster = cluster_labels_12
) %>%
  arrange(cluster, taxon)

cluster_df_12

write.csv(cluster_df_12, file.path(microbiome_dir, "GMCcol_taxa_cluster_assignments_k12.csv"), row.names = FALSE)



#Final heatmap
# --- Order taxa by manual cluster assignment ---

ordered_taxa <- names(sort(cluster_labels_12))
cor_mat_ordered <- cor_mat_all[ordered_taxa, ordered_taxa]

# --- Annotation bar showing cluster membership ---
annotation_df <- data.frame(
  Cluster = factor(cluster_labels_12[ordered_taxa])
)
rownames(annotation_df) <- ordered_taxa

# --- Compute gap positions: where cluster ID changes along the ordered taxa ---
cluster_order <- cluster_labels_12[ordered_taxa]
gaps <- which(diff(as.numeric(factor(cluster_order))) != 0)

# --- Color palette (reuse your existing one) ---
cor_cols <- colorRampPalette(c("blue", "white", "red"))(100)
cor_breaks <- seq(-1, 1, length.out = 101)

pdf(
  file.path(microbiome_dir, "taxa_heatmap_12_clusters.pdf"),
  width = 22,
  height = 20
)
pheatmap(
  cor_mat_ordered,
  color = cor_cols,
  breaks = cor_breaks,
  cluster_rows = FALSE,   # turn off auto-clustering; we're using manual order
  cluster_cols = FALSE,
  gaps_row = gaps,
  gaps_col = gaps,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  fontsize_row = 12,
  fontsize_col = 12,
  main = "Taxa Spearman correlation heatmap | Cluster assignment"
)
dev.off()


# Make covariate table with all 136 taxa

# --- make covariate file ---
rep_taxa_clr <- clr_sub %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleName")
  
rep_taxa_final <- colnames(rep_taxa_clr)

colnames(rep_taxa_clr) <- c("SampleName", gsub("[ /]", "_", rep_taxa_final))

cv_with_taxa <- cv %>%
  filter(Treatment == "microbiome") %>%
  left_join(rep_taxa_clr, by = "SampleName")

write.table(cv_with_taxa,
            file.path(cov_dir, "GMCcol_deep_covariates_all_microbial_taxa_7_29_2026.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

