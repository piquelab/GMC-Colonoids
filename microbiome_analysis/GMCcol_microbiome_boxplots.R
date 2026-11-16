library(tidyverse)
library(ggplot2)

microbiome_dir <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/microbiome_analysis"
cov.file <- "/wsu/home/groups/piquelab/colonoids/GMCcol/GMCcol_deep/covariates/GMCcol_deep_covariates.txt"

cv <- read.table(cov.file, sep="\t", stringsAsFactors=FALSE, header=TRUE, comment="")

cv_microbiome <- cv %>%
filter(Urbanism %in% c("Urban","Rural")) %>%
distinct(Microbiome, Urbanism, Prevotella, Bacteroides, Firmicutes) %>%
mutate(microbiome_type=ifelse(Urbanism == "Urban", "Urban Microbiome", "Rural Microbiome"),
microbiome_type=factor(microbiome_type, levels=c("Urban Microbiome","Rural Microbiome")))

mycolors <- c("Urban Microbiome"="#0072B2", "Rural Microbiome"="#009E73")

p_label <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    "ns"
}

plot_microbiome_abundance <- function(data, taxon) {

    tmp <- data %>%
    select(microbiome_type, abundance=all_of(taxon)) %>%
    drop_na()
    
    urban <- tmp$abundance[tmp$microbiome_type == "Urban Microbiome"]
    rural <- tmp$abundance[tmp$microbiome_type == "Rural Microbiome"]
    
    alternative <- ifelse(taxon == "Prevotella", "less", "greater")
    test_res <- wilcox.test(urban, rural, alternative=alternative, exact=TRUE)
    
    pval <- test_res$p.value
    sig <- p_label(pval)
    pval_text <- ifelse(pval < 0.001, "P < 0.001", paste0("P = ", formatC(pval, format="f", digits=3)))
    
    y_max <- max(tmp$abundance, na.rm=TRUE)
    y_min <- min(tmp$abundance, na.rm=TRUE)
    y_range <- y_max-y_min
    if (y_range == 0) y_range <- 1
    
    y_bracket <- y_max+0.12*y_range
    tick <- 0.04*y_range
    
    p <- ggplot(tmp, aes(x=microbiome_type, y=abundance, color=microbiome_type)) +
    geom_boxplot(width=0.5, fill=NA, outlier.shape=NA, linewidth=0.8) +
    geom_point(size=4, alpha=0.9) +
    annotate("segment", x=1, xend=2, y=y_bracket, yend=y_bracket, linewidth=0.7) +
    annotate("segment", x=1, xend=1, y=y_bracket, yend=y_bracket-tick, linewidth=0.7) +
    annotate("segment", x=2, xend=2, y=y_bracket, yend=y_bracket-tick, linewidth=0.7) +
    annotate("text", x=1.5, y=y_bracket+0.05*y_range, label=sig, size=6) +
    annotate("text", x=1.5, y=y_bracket+0.13*y_range, label=pval_text, size=5) +
    scale_color_manual(values=mycolors) +
    scale_x_discrete(labels=c("Urban Microbiome"="Urban\nMicrobiome", "Rural Microbiome"="Rural\nMicrobiome")) +
    scale_y_continuous(expand=expansion(mult=c(0.08,0.30))) +
    labs(x=NULL, y=bquote(italic(.(taxon))~"abundance")) +
    theme_classic(base_size=18) +
    theme(axis.text.x=element_text(size=16, color="black", hjust=0.5, lineheight=0.9),
    axis.text.y=element_text(size=16, color="black"),
    axis.title=element_text(size=18),
    legend.position="none")
    
    safe_taxon <- gsub("[^A-Za-z0-9._-]", "_", taxon)
    
    ggsave(file.path(microbiome_dir, paste0("GMCcol_microbiome_samples_", safe_taxon, "_boxplot.png")), p, width=5, height=5, dpi=600, bg="white")
    ggsave(file.path(microbiome_dir, paste0("GMCcol_microbiome_samples_", safe_taxon, "_boxplot.pdf")), p, width=5, height=5, bg="white")
    
    cat(taxon, "- Wilcoxon rank-sum P =", pval, "\n")
    
    return(p)
}


p_prevotella <- plot_microbiome_abundance(cv_microbiome, "Prevotella")
p_bacteroides <- plot_microbiome_abundance(cv_microbiome, "Bacteroides")
p_firmicutes <- plot_microbiome_abundance(cv_microbiome, "Firmicutes")
