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

# Categorise HDL levels (mmol/L)
meta$HDL_cat <- cut(meta$HDL,
                    breaks = c(-Inf, 1.0, 1.5, Inf),
                    labels = c("Low", "Desirable", "Optimal"),
                    right = FALSE)
meta$HDL_cat <- factor(meta$HDL_cat, levels = c("Low", "Desirable", "Optimal"))

# Compute Bray-Curtis distance
bray_dist <- vegdist(otu_mat, method = "bray")

# PERMANOVA
permanova_hdl <- adonis2(bray_dist ~ HDL_cat, data = meta)
print(permanova_hdl)

# Extract F and p-values
permanova_F <- round(permanova_hdl$F[1], 3)
permanova_p <- permanova_hdl$`Pr(>F)`[1]

# PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
pcoa_df <- data.frame(pcoa$points, HDL_cat = meta$HDL_cat)
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

plot_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = HDL_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = paste0("PCoA (Bray-Curtis)\nPERMANOVA F=", permanova_F, ", p=", permanova_p),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  )

# PCA
pca <- prcomp(otu_mat, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], HDL_cat = meta$HDL_cat)
var_pca <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = HDL_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA by HDL category",
    x = paste0("PC1 (", var_pca[1], "%)"),
    y = paste0("PC2 (", var_pca[2], "%)")
  )

# NMDS
nmds <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(nmds$points, HDL_cat = meta$HDL_cat)

plot_nmds <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = HDL_cat)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis) by HDL category")

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
     tip.color = as.numeric(meta$HDL_cat),
     main = "UPGMA Tree (Bray-Curtis) by HDL category")

# ANOSIM
anosim_hdl <- anosim(bray_dist, grouping = meta$HDL_cat)
print(anosim_hdl)

# ANOSIM R and p-value
anosim_R <- round(anosim_hdl$statistic, 3)
anosim_p <- anosim_hdl$signif

# Distance data for ggplot boxplot
anosim_dist <- as.matrix(bray_dist)
group_labels <- meta$HDL_cat

get_group_type <- function(i, j) {
  if (group_labels[i] == group_labels[j]) {
    return("Within")
  } else {
    return("Between")
  }
}

group_comparison <- combn(seq_along(group_labels), 2, function(idx) {
  i <- idx[1]
  j <- idx[2]
  data.frame(
    Distance = anosim_dist[i, j],
    Group = get_group_type(i, j)
  )
}, simplify = FALSE)

anosim_df <- bind_rows(group_comparison)

# ANOSIM Boxplot
ggplot(anosim_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste0("ANOSIM by HDL category\nR = ", anosim_R, ", p = ", anosim_p))
