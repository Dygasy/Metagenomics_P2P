# === Load libraries ===
library(readxl)
library(tidyverse)
library(pheatmap)

# === Load metadata and sort by BMI ===
metadata <- read_excel("E:/Krona_results/05_Taxonomy/Metadata/Clinical/Meta_Data_Stool_2024-12-06.xlsm")
metadata_sorted <- metadata %>%
  arrange(desc(BMI))  # Replace with correct BMI column name

# === Define taxonomic levels and file paths ===
levels <- c("Phylum", "Class", "Order", "family", "Genus", "Species")
base_path <- "E:/Krona_results/05_Taxonomy/Relative abundance"
file_paths <- paste0(base_path, "/", levels, "_relativeabundance.txt")

# === Loop through each taxonomic level ===
for (i in seq_along(levels)) {
  level <- levels[i]
  file <- file_paths[i]
  
  # Load abundance data
  df <- read_tsv(file, show_col_types = FALSE)
  colnames(df)[1] <- "Taxon"
  
  # === Select top 20 most abundant taxa ===
  top_taxa <- df %>%
    pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
    group_by(Taxon) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
    slice_max(Total, n = 20) %>%
    pull(Taxon)
  
  df_top <- df %>% filter(Taxon %in% top_taxa)
  
  # === Transpose for heatmap ===
  df_long <- df_top %>%
    pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance") %>%
    pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)
  
  # === Merge with metadata ===
  merged <- metadata_sorted %>%
    select(SampleID = `Participant_ID`, BMI) %>%  # Replace with actual sample ID column
    inner_join(df_long, by = "SampleID")
  
  # === Prepare matrix ===
  rownames_mat <- merged$SampleID
  heatmap_data <- merged %>% select(-SampleID, -BMI)
  rownames(heatmap_data) <- rownames_mat
  
  # === Plot heatmap ===
  pheatmap(
    mat = as.matrix(heatmap_data),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = paste("GTDB", level, "- Top 20 Abundance (Sorted by BMI)"),
    fontsize = 8,
    scale = "row"
  )
}
