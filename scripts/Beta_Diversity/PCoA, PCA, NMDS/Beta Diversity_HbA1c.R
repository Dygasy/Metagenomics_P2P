# Load libraries
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Extract OTU table and metadata
otu_mat <- as(otu_table(ps), "matrix")
otu_mat <- t(otu_mat)  # Samples as rows
meta <- as(sample_data(ps), "data.frame")

# Categorize HbA1c
meta$HbA1c_cat <- cut(
  meta$HbA1c,
  breaks = c(-Inf, 5.9, 6.4, Inf),
  labels = c("Normal", "Prediabetes", "Diabetes")
)

# Ensure it's a factor
meta$HbA1c_cat <- factor(meta$HbA1c_cat, levels = c("Normal", "Prediabetes", "Diabetes"))

# Compute Bray-Curtis distance
bray_dist <- vegdist(otu_mat, method = "bray")

# PERMANOVA
permanova_hba1c <- adonis2(bray_dist ~ HbA1c_cat, data = meta)
print(permanova_hba1c)
permanova_F <- round(permanova_hba1c$F[1], 3)
permanova_p <- permanova_hba1c$`Pr(>F)`[1]

# PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
pcoa_df <- data.frame(pcoa$points, HbA1c_cat = meta$HbA1c_cat)
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

plot_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HbA1c_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = paste0("PCoA (Bray-Curtis)\nPERMANOVA F=", permanova_F, ", p=", permanova_p),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  )

# PCA
pca <- prcomp(otu_mat, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], HbA1c_cat = meta$HbA1c_cat)
var_pca <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = HbA1c_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA by HbA1c Category",
    x = paste0("PC1 (", var_pca[1], "%)"),
    y = paste0("PC2 (", var_pca[2], "%)")
  )

# NMDS
nmds <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(nmds$points, HbA1c_cat = meta$HbA1c_cat)

plot_nmds <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = HbA1c_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis) by HbA1c Category")

# Combine ordination plots
ggarrange(
  plot_pcoa, plot_pca, plot_nmds,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "right"
)

# UPGMA Tree
tree <- hclust(bray_dist, method = "average")
plot(as.phylo(tree), type = "fan",
     tip.color = as.numeric(meta$HbA1c_cat),
     main = "UPGMA Tree (Bray-Curtis) by HbA1c Category")

# ANOSIM
anosim_result <- anosim(bray_dist, grouping = meta$HbA1c_cat)
print(anosim_result)
anosim_R <- round(anosim_result$statistic, 3)
anosim_p <- anosim_result$signif

# ANOSIM boxplot
anosim_dist <- as.matrix(bray_dist)
group_labels <- meta$HbA1c_cat

get_group_type <- function(i, j) {
  if (group_labels[i] == group_labels[j]) {
    return("Within")
  } else {
    return("Between")
  }
}

group_comparison <- combn(seq_along(group_labels), 2, function(idx) {
  i <- idx[1]; j <- idx[2]
  data.frame(Distance = anosim_dist[i, j], Group = get_group_type(i, j))
}, simplify = FALSE)

anosim_df <- bind_rows(group_comparison)

ggplot(anosim_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste0("ANOSIM by HbA1c Category\nR = ", anosim_R, ", p = ", anosim_p))
