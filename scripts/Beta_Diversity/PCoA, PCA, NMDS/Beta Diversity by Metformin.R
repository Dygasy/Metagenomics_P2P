# -------------------------------
# Load libraries
library(phyloseq)
library(vegan)
library(ape)
library(ggplot2)
library(dplyr)
library(ggpubr)

# -------------------------------
# Extract OTU table and metadata
otu_mat <- as(otu_table(ps), "matrix")
otu_mat <- t(otu_mat)  # Samples as rows
meta <- as(sample_data(ps), "data.frame")

# -------------------------------
# Ensure Metformin is a factor
meta$Metformin <- factor(meta$Metformin, levels = c("no", "yes"))

# -------------------------------
# Compute Bray-Curtis distance
bray_dist <- vegdist(otu_mat, method = "bray")

# -------------------------------
# PERMANOVA
permanova <- adonis2(bray_dist ~ Metformin, data = meta)
print(permanova)

# Extract PERMANOVA F and p-values
permanova_F <- round(permanova$F[1], 3)
permanova_p <- permanova$`Pr(>F)`[1]

# -------------------------------
# PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
pcoa_df <- data.frame(pcoa$points, Metformin = meta$Metformin)
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

plot_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Metformin)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = paste0("PCoA (Bray-Curtis) by Metformin\nPERMANOVA F=", permanova_F, ", p=", permanova_p),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  )

# -------------------------------
# PCA
pca <- prcomp(otu_mat, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], Metformin = meta$Metformin)
var_pca <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Metformin)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA by Metformin",
    x = paste0("PC1 (", var_pca[1], "%)"),
    y = paste0("PC2 (", var_pca[2], "%)")
  )

# -------------------------------
# NMDS
nmds <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(nmds$points, Metformin = meta$Metformin)

plot_nmds <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Metformin)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis) by Metformin")

# -------------------------------
# Combine ordination plots
ggarrange(
  plot_pcoa, plot_pca, plot_nmds,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "right"
)

# -------------------------------
# UPGMA Tree
tree <- hclust(bray_dist, method = "average")
plot(as.phylo(tree), type = "fan",
     tip.color = as.numeric(meta$Metformin),
     main = "UPGMA Tree (Bray-Curtis) by Metformin")

# -------------------------------
# ANOSIM
anosim_res <- anosim(bray_dist, grouping = meta$Metformin)
print(anosim_res)

# Extract ANOSIM R and p-value
anosim_R <- round(anosim_res$statistic, 3)
anosim_p <- anosim_res$signif

# -------------------------------
# Construct data frame for ggplot
anosim_dist <- as.matrix(bray_dist)
group_labels <- meta$Metformin

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

# -------------------------------
# ANOSIM plot
ggplot(anosim_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste0("ANOSIM by Metformin\nR = ", anosim_R, ", p = ", anosim_p))
