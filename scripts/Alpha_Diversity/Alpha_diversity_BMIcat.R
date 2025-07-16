# -------------------------------
# Load necessary libraries
library(tidyverse)
library(readxl)
library(phyloseq)
library(vegan)
library(ggpubr)

# -------------------------------
# 1. Load OTU table
otu_raw <- read.delim("E:/Krona_results/05_Taxonomy/Absolute Abundance/Species_abundance.txt",
                      header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Clean sample IDs in column names
colnames(otu_raw) <- gsub("_diamond_output_.*", "", colnames(otu_raw))
otu_mat <- as.matrix(otu_raw)

# -------------------------------
# 2. Load metadata
metadata <- read_excel("E:/Krona_results/Metadata/Clinical/Meta_Data_Stool_2024-12-06.xlsm",
                       sheet = "data") %>% 
  as.data.frame()

# Set rownames as Participant ID
rownames(metadata) <- metadata$`Participant_ID`

# -------------------------------
# 3. Match sample IDs between metadata and OTU table
shared_ids <- intersect(colnames(otu_mat), rownames(metadata))
otu_mat <- otu_mat[, shared_ids]
metadata <- metadata[shared_ids, ]

# -------------------------------
# 4. Create phyloseq object
otu_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)
sample_ps <- sample_data(metadata)
ps <- phyloseq(otu_ps, sample_ps)

# -------------------------------
# 5. Estimate alpha diversity
alpha_div <- estimate_richness(ps, measures = c("Chao1", "Shannon", "Simpson"))

# Add Fisher alpha diversity manually
otu_counts <- t(as(otu_table(ps), "matrix"))
fisher_values <- apply(otu_counts, 1, fisher.alpha)
alpha_div$Fisher <- fisher_values

# Merge with metadata
alpha_div$SampleID <- rownames(alpha_div)
metadata$SampleID <- rownames(metadata)
alpha_combined <- left_join(alpha_div, metadata, by = "SampleID")

# -------------------------------
# 6. Fix NA and add Underweight label
alpha_combined$BMIcat[is.na(alpha_combined$BMIcat)] <- "Underweight"

# Reorder factor levels
alpha_combined$BMIcat <- factor(alpha_combined$BMIcat,
                                levels = c("Underweight",
                                           "Normal - Low Risk (healthy range)",
                                           "Overweight - Moderate Risk",
                                           "Obese - High Risk"))

# -------------------------------
# 7. Set pairwise comparison groups
pairwise_comparisons <- list(
  c("Underweight", "Normal - Low Risk (healthy range)"),
  c("Underweight", "Overweight - Moderate Risk"),
  c("Underweight", "Obese - High Risk"),
  c("Normal - Low Risk (healthy range)", "Overweight - Moderate Risk"),
  c("Normal - Low Risk (healthy range)", "Obese - High Risk"),
  c("Overweight - Moderate Risk", "Obese - High Risk")
)

# -------------------------------
# 8. Plotting function for diversity metrics
plot_metric <- function(data, metric, y_label) {
  ggplot(data, aes(x = BMIcat, y = .data[[metric]], fill = BMIcat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity"), y = y_label, x = "BMI Category") +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format")
}

# Generate all plots
plot_chao1   <- plot_metric(alpha_combined, "Chao1", "Chao1")
plot_shannon <- plot_metric(alpha_combined, "Shannon", "Shannon")
plot_simpson <- plot_metric(alpha_combined, "Simpson", "Simpson")
plot_fisher  <- plot_metric(alpha_combined, "Fisher", "Fisher")

# -------------------------------
# 9. Combine and display all plots
ggarrange(plot_chao1, plot_shannon, plot_simpson, plot_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
