Functional pathway to clinical metadata:
1. aggregate this KEGG pathway counts matrix per sample.
2. use differential testing approach (MaAslin3, ALDEx2, DESeq2 or ANOVA) on the pathway abundance data across your metadata categories to identify differentially enriched pathways across categories.
3. this is to link to metadata for testing significance across categories

For now, we have the following data: 
1. KEGG_pathway_counts_per_sample.tsv
2. metadata_lipid_numeric.tsv
3. merge them into one file before modelling or visualisation

```bash
# Load required libraries
library(dplyr)
library(readr)
library(tibble)

# Step 1: Load pathway counts and metadata
kegg_counts <- read_tsv("KEGG_pathway_counts_per_sample.tsv") %>%
  column_to_rownames("Sample")

metadata_tmp$Sample <- rownames(kegg_counts)
metadata <- metadata_tmp %>% column_to_rownames("Sample")

# Step 2: Match and combine based on sample IDs
common_samples <- intersect(rownames(kegg_counts), rownames(metadata))
kegg_counts <- kegg_counts[common_samples, ]
metadata <- metadata[common_samples, ]

# Step 3: Merge metadata and pathway data
combined_df <- cbind(metadata, kegg_counts)

# Step 4: Correlation Analysis
cor_matrix <- cor(combined_df[, c("LDL", "HDL", "Cholesterol")], 
                  combined_df[, grep("ko", colnames(combined_df))], 
                  use = "pairwise.complete.obs")
pheatmap(cor_matrix, main = "Correlation: Lipids vs KEGG Pathways")

# Step 5: Export for MaAsLin3 or visualization
write.table(merged_df, "merged_kegg_pathway_metadata.tsv", sep = "\t", quote = FALSE)
```
Step 6: Linear Model Example - LDL vs Top Pathways
```bash
lm_results <- lapply(colnames(kegg_counts), function(pathway) {
  fit <- lm(combined_df[[pathway]] ~ combined_df$LDL)
  summary(fit)$coefficients[2, ]  # Extract slope info for LDL
})
#Convert matrix to a data frame explicitly
lm_df <- do.call(rbind, lm_results)
rownames(lm_df) <- colnames(kegg_counts)
colnames(lm_df) <- c("Estimate", "Std.Error", "t.value", "p.value")

# View significant results
significant_kegg <- lm_df[lm_df[, "p.value"] < 0.05, ]
head(significant_kegg[order(significant_kegg[, "p.value"]), ])

#Step 7: Linear Models for all three lipid traits 
traits <- c("LDL", "HDL", "Cholesterol")
lm_all <- list()

for (trait in traits) {
  lm_results <- lapply(colnames(kegg_counts), function(pathway) {
    fit <- lm(combined_df[[pathway]] ~ combined_df[[trait]])
    summary(fit)$coefficients[2, ]
  })
  
  lm_df <- as.data.frame(do.call(rbind, lm_results))  # Convert to data frame
  rownames(lm_df) <- colnames(kegg_counts)
  colnames(lm_df) <- c("Estimate", "Std.Error", "t.value", "p.value")
  
  lm_df$Trait <- trait
  lm_df$Pathway <- rownames(lm_df)
  
  lm_all[[trait]] <- lm_df
}

lm_all_df <- do.call(rbind, lm_all)
# Step 8: Filter significant results (p<0.05)
sig_results <- lm_all_df %>% 
  filter(p.value < 0.05) %>% 
  arrange(p.value)

# Step 9: Barplot top associations
library(ggplot2)

ggplot(sig_results, aes(x = reorder(Pathway, Estimate), y = Estimate, fill = Trait)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~Trait, scales = "free_y") +
  labs(title = "Significant KEGG Pathway Associations with Lipid Traits",
       x = "KEGG Pathway", y = "Effect Size (Estimate)") +
  theme_minimal()
Step 10: Heatmap of Coefficients
library(pheatmap)

# Spread into wide format for heatmap
coef_mat <- sig_results %>%
  select(Pathway, Trait, Estimate) %>%
  pivot_wider(names_from = Trait, values_from = Estimate, values_fill = 0) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

pheatmap(coef_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap: KEGG Pathways vs Lipid Traits")
# interpretation
✅ Heatmap: Correlation of KEGG Pathways vs Lipids
Each row represents a KEGG pathway (ko ID).

Each column is a lipid trait (LDL, HDL, Cholesterol).

Color gradient: Red indicates a stronger positive association (higher β estimate), blue indicates negative association.

From the plot, most correlations are small and negative (blue), except for a few cholesterol-related pathways that stand out slightly (red).

✅ Barplot: Significant KEGG Pathway Associations
Filters only those with p < 0.05.

The top hit seems to be ko04090 and ko04052, both associated with Cholesterol.

The bar height reflects the effect size (Estimate from linear model).
