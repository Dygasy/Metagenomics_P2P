GTDB to produce the following outputs after MEGAN processing: 
1. Relative abundances from different taxonomic levels (Phylum, Class, Order, Genus, Species). 
2. Absolute abundances from different taxonomic levels (Phylum, Class, Order, Genus, Species),

Coupled with the extensive list of clinical metadata that's being provided. 
install the following packages in R
```bash
# Core data handling & visualization
install.packages(c("tidyverse", "readxl", "ggpubr"))

# Microbiome-specific
install.packages(c("phyloseq", "microbiome", "vegan"))

# For advanced plots
install.packages(c("ggplot2", "pheatmap", "ComplexHeatmap", "ggrepel"))

# install_packages.R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "DESeq2", "edgeR"))
install.packages(c("vegan", "microbiome", "ggplot2", "readxl", "pheatmap", "ComplexHeatmap"))

```

changed my otu files into a csv file
*seek out a pre-define small set of primary clinical variables of interest, based on biological hypotheses or clinical relevance. Sometimes secondary factors can also be considered. 
generate box plots + p-values for categorical and correlation for continuous. 
If many variables, these can go into Supplementary Figures. 

NHCS automated alpha diversity vs metadata
2-5 clinical or exposure factors, boxplots + p-values (Wilcoxon, Kruskal-Wallis)
```bash
# ==============================================
#  NHCS automated alpha diversity vs metadata
# ==============================================

library(tidyverse)
library(ggpubr)

# Assuming `alpha_combined` exists (alpha diversity + metadata joined)
alpha_metrics <- c("Chao1", "Shannon", "Simpson", "Fisher")

# Use a smaller list or your selected important variables
clinical_vars <- c("BMI", "SBP", "DBP", "Gender", "Weight", "Height",
                   "Pulse", "WaistCircumference", "Hypertension",
                   "Dyslipidemia", "Diabetes_mellitus", "Smoking",
                   "ARB", "ACE_inhibitor", "Calcium_channel_blocker",
                   "Alcohol", "Age", "VO2Max")

# Helper function for each plot
multi_alpha_plot <- function(var) {
  plots <- list()
  
  for (metric in alpha_metrics) {
    if (var %in% colnames(alpha_combined)) {
      p <- NULL
      
      if (is.numeric(alpha_combined[[var]])) {
        # Continuous
        p <- ggscatter(
          alpha_combined, x = var, y = metric,
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          title = paste(metric, "vs", var)
        ) + theme_bw()
        
      } else {
        # Categorical
        unique_vals <- unique(na.omit(alpha_combined[[var]]))
        test_method <- if(length(unique_vals) == 2) "wilcox.test" else "kruskal.test"
        p <- ggplot(alpha_combined, aes(x = .data[[var]], y = .data[[metric]], fill = .data[[var]])) +
          geom_boxplot() + geom_jitter(width=0.2, alpha=0.6) +
          stat_compare_means(method = test_method, label = "p.format") +
          theme_bw() +
          labs(title = paste(metric, "by", var), x=var, y=metric)
      }
      
      plots[[metric]] <- p
    }
  }
  
  return(plots)
}

# Directory to save
dir.create("alpha_multi_metric_plots", showWarnings = FALSE)

# Loop through all clinical vars
for (var in clinical_vars) {
  plots_list <- multi_alpha_plot(var)
  # If there are plots, combine them into a 2x2 grid
  if (length(plots_list) > 0) {
    combined_plot <- ggarrange(plotlist = plots_list, ncol=2, nrow=2, common.legend = TRUE)
    ggsave(filename = paste0("alpha_multi_metric_plots/alpha_vs_", var, ".png"),
           plot = combined_plot, width = 10, height = 8, dpi = 300)
  }
}

```

Perform Beta Diversity (ordination) and PERMANOVA for overall microbiome composition differ by the main grouping factors. 
you can adjust for covariates (via marginal or multivariate PERMANOVA)
Based on your metadata variables, spot the one with the most biologically or clinically relevant to show in your main figures. The rest can go to the Supplementary Figures
or an exploratory data appendix.
If there are many comparisons to be made, you will have to adjust p values (eg; False Discovery Rate , Bonferroni)
Theres a difference between exploratory alpha diversity and confirmatory. Exploratory alpha diversity checks are usually carefully framed. Eg Ëxploratory correlations were performed
between diversity metrics and continuous clinical variables, with p-values adjusted by Benjamini-Hochberg.

Include: (PERMANOVA+ PCoA ordinations + Adonis tests) 
if you test on every clinical parameter, it will create false positives (FDR) because of multiple testing burden and clutters the story where its hard to interpret dozens of comparisons. 

PERMANOVA: standard non parametric way to test whether microbial communities differ in composition across groups or along gradients
* no assumption for multivariate normality
* use permutation to generate p-values
* computes an r^2 which is an interpretable effect size to show how much variation in community structure is explained by your factor.
R-script to automate Beta diversity (PERMANOVA with Adonis2)
```bash
# ===========================================
# 🚀 Automated one-by-one PERMANOVA for beta diversity
# ===========================================

```bash
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)

# Distance matrix
bray_dist <- phyloseq::distance(physeq, method = "bray")
ordination <- ordinate(physeq, method = "PCoA", distance = bray_dist)
ordination_df <- as.data.frame(ordination$vectors[,1:2])
ordination_df$SampleID <- rownames(ordination_df)
ordination_df <- cbind(ordination_df, meta_df[rownames(ordination_df), ])

# Factors to plot
factors_to_plot <- c("Gender", "Hypertension", "BMI", "Age")

# Loop over factors
for (factor in factors_to_plot) {
  
  # PERMANOVA
  keep_samples <- complete.cases(meta_df[[factor]])
  this_dist <- as.dist(as.matrix(bray_dist)[keep_samples, keep_samples])
  this_meta <- meta_df[keep_samples,, drop=FALSE]
  formula <- as.formula(paste("this_dist ~", factor))
  adonis_out <- adonis2(formula, data = this_meta, permutations = 999)
  p_val <- adonis_out$`Pr(>F)`[1]
  R2 <- adonis_out$R2[1]
  
  # Prepare plot data
  plot_df <- ordination_df[keep_samples, ]
  
  # Plot
  p <- ggplot(plot_df, aes_string(x="Axis.1", y="Axis.2", color=factor)) +
    geom_point(size=3, alpha=0.8) +
    stat_ellipse(type="norm", linetype=2) +
    theme_minimal() +
    labs(
      title=paste0("PCoA (Bray-Curtis) colored by ", factor),
      subtitle=sprintf("PERMANOVA p=%.3f, R²=%.3f", p_val, R2),
      x="PCoA1", y="PCoA2", color=factor
    ) +
    theme(legend.position = "right")
  
  print(p)
}

```
