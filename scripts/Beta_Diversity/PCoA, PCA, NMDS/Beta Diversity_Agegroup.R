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

# Create Age_group column based on Age
meta$Age_group <- cut(meta$Age,
                      breaks = c(-Inf, 39, 59, Inf),
                      labels = c("18–39", "40–59", "60 and above"),
                      right = TRUE)

# Make sure it's a factor with correct order
meta$Age_group <- factor(meta$Age_group, levels = c("18–39", "40–59", "60 and above"))

# Compute Bray-Curtis distance
bray_dist <- vegdist(otu_mat, method = "bray")

# PERMANOVA by Age_group
permanova_age <- adonis2(bray_dist ~ Age_group, data = meta)
print(permanova_age)

# Extract F and p-values
permanova_F <- round(permanova_age$F[1], 3)
permanova_p <- permanova_age$`Pr(>F)`[1]

# PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 1)
pcoa_df <- data.frame(pcoa$points, Age_group = meta$Age_group)
colnames(pcoa_df)[1:2] <- c("PCoA1", "PCoA2")

plot_pcoa <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Age_group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = paste0("PCoA (Bray-Curtis)\nPERMANOVA F=", permanova_F, ", p=", permanova_p),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  )

# PCA
pca <- prcomp(otu_mat, scale. = TRUE)
pca_df <- data.frame(pca$x[, 1:2], Age_group = meta$Age_group)
var_pca <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)

plot_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Age_group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA by Age Group",
    x = paste0("PC1 (", var_pca[1], "%)"),
    y = paste0("PC2 (", var_pca[2], "%)")
  )

# NMDS
nmds <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(nmds$points, Age_group = meta$Age_group)

plot_nmds <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Age_group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS (Bray-Curtis) by Age Group")

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
     tip.color = as.numeric(meta$Age_group),
     main = "UPGMA Tree (Bray-Curtis) by Age Group")

# ANOSIM
anosim_age <- anosim(bray_dist, grouping = meta$Age_group)
print(anosim_age)

# ANOSIM R and p-value
anosim_R <- round(anosim_age$statistic, 3)
anosim_p <- anosim_age$signif

# Construct data frame for ggplot ANOSIM boxplot
anosim_dist <- as.matrix(bray_dist)
group_labels <- meta$Age_group

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

# ANOSIM ggplot boxplot
ggplot(anosim_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = paste0("ANOSIM by Age Group\nR = ", anosim_R, ", p = ", anosim_p))
