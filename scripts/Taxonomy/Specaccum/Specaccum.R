library(tidyverse)
library(vegan)

# === 1. Load your species abundance (absolute count) table ===
df <- read_tsv("E:/Krona_results/05_Taxonomy/Absolute Abundance/Species_abundance.txt")

# === 2. Convert to matrix with species as rows ===
df_matrix <- df %>%
  column_to_rownames(var = colnames(df)[1]) %>%
  as.matrix()

# === 3. Transpose: samples as rows ===
df_matrix <- t(df_matrix)

# === 4. Generate species accumulation object ===
spec_accum <- specaccum(df_matrix, method = "random")

# === 5. Extract species richness matrix from permutations ===
spec_matrix <- spec_accum$perm

# === 6. Prepare data frame for ggplot ===
spec_df <- as.data.frame(spec_matrix)
spec_df$SampleOrder <- 1:nrow(spec_df)
spec_long <- pivot_longer(spec_df, cols = -SampleOrder, names_to = "Permutation", values_to = "Species")

# === 7. Plot as boxplot with ggplot ===
ggplot(spec_long, aes(x = SampleOrder, y = Species)) +
  geom_boxplot(fill = "steelblue", outlier.size = 0.5, alpha = 0.7) +
  labs(
    title = "Species Accumulation Curve (with boxplots)",
    x = "Number of sample",
    y = "Observed species"
  ) +
  theme_minimal(base_size = 14)
