# === Load required libraries ===
library(tidyverse)
library(readxl)
library(pheatmap)

# === 1. Load metadata ===
metadata <- read_excel("E:/Krona_results/05_Taxonomy/Metadata/Clinical/Meta_Data_Stool_2024-12-06.xlsm")

metadata <- metadata %>%
  rename(SampleID = Participant_ID) %>%
  select(SampleID, BMI) %>%
  drop_na(BMI)

# === 2. Load Phylum-level relative abundance ===
abundance <- read_tsv("E:/Krona_results/05_Taxonomy/Relative abundance/Phylum_relativeabundance.txt", show_col_types = FALSE)
colnames(abundance)[1] <- "Taxon"

# === 3. Identify top 20 most abundant genera ===
top_genera <- abundance %>%
  pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Taxon) %>%
  summarise(total = sum(Abundance, na.rm = TRUE)) %>%
  slice_max(total, n = 20) %>%
  pull(Taxon)

# === 4. Filter for top genera only ===
abundance_top <- abundance %>% filter(Taxon %in% top_genera)

# === 5. Reshape and merge with BMI ===
abund_matrix <- abundance_top %>%
  pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance") %>%
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)

merged <- metadata %>%
  inner_join(abund_matrix, by = "SampleID")

# === 6. Prepare numeric matrix ===
heatmap_data <- merged %>%
  select(-SampleID, -BMI)

# 🔒 Force numeric conversion, then check
heatmap_data <- heatmap_data %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  as.data.frame()

# Print to confirm all numeric
print(sapply(heatmap_data, class))  # All must be "numeric"

# Set rownames
rownames(heatmap_data) <- merged$SampleID

# === 7. Prepare BMI annotation ===
annotation_row <- merged %>% select(BMI) %>% as.data.frame()
rownames(annotation_row) <- merged$SampleID

# === 8. Generate heatmap ===
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_row = annotation_row,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  main = "Top 20 PhylumAbundance (Clustered by Microbiome Profile, Annotated by BMI)",
  scale = "row",
  fontsize_row = 7,
  fontsize_col = 8
)

