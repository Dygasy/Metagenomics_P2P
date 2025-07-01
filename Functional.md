Run emapper.py (eggNOG-mapper v2)
1. Make sure your environment is activated
```bash
conda activate eggnogv2
```
create directory for eggNOG database properly
```bash
mkdir -p /home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data
```
```bash
download_eggnog_data.py --data_dir /home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data
```
uncompress the eggnog_proteins.dmnd.gz before emapper.py can be used.

Use of .faa file from gmhmmmp_*.faa file to run eggNOG-mapper for each sample! batch script
```bash
nano run_eggnog_mapper.sh
```
```bash
#!/bin/bash

# Base directory containing AA_* folders
base_dir="/mnt/e/Krona_results/below"

# Output directory for EggNOG results
output_dir="/mnt/e/Krona_results/eggnog_mapper_output"
mkdir -p "$output_dir"

# Path to EggNOG-mapper data directory (update if yours is different)
data_dir="/home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data"


# Loop through all subfolders with FAA files
for folder in "$base_dir"/AA_*; do
    faa_file=$(find "$folder" -type f -name "*.faa")
    
    if [ -f "$faa_file" ]; then
        folder_name=$(basename "$folder")
        output_prefix="${output_dir}/${folder_name}_eggnog"

        echo "✅ Processing: $faa_file"
        
        emapper.py \
            -i "$faa_file" \
            --itype proteins \
            --output "$output_prefix" \
            --output_dir "$output_dir" \
            --data_dir "$data_dir" \
            --cpu 4
    else
        echo "❌ No FAA file found in $folder"
    fi
done
```
Using DIAMOND v2.1.6 to run emapper.py, run this script. 
output -- *_eggnog.emapper.annotations, *_eggnog.emapper.hits, *_eggnog.emapper.seed_orthologs
eggnog emapper annotations provides main functional annotation output. to be used for pathway-level analysis, GO terms, KEGG, eggNOG categories (COGs), etc
emapper.hits provides raw DIAMOND alignment (can be ignored unless troubleshooting)
emapper_seed_orthologs provide closest ortholog matches - to be used for fine-level comparative analysis

For the functional analysis, we would like to focus on:
1. KEGG orthology (KO)
2. eggNOG COG categories
3. Focus on GO terms

Extract KO Terms from .annotations files
R script:
```bash
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

path <- "E:/eggnog_mapper_output/"
files <- list.files(path, pattern = "\\.annotations$", full.names = TRUE)

kegg_list <- list()

for (f in files) {
  sample_id <- str_remove(basename(f), "_eggnog.emapper.annotations")
  
  df <- read_tsv(f, skip = 4, comment = "#", show_col_types = FALSE)
  colnames(df) <- make.names(colnames(df), unique = TRUE)
  
  # Extract KO and KEGG Pathway columns (column 12 = KO, column 17 = Pathway based on your screenshot)
  ko_present <- grepl("^ko:", df[[12]])
  path_present <- grepl("^ko", df[[17]])
  
  ko_terms <- df[ko_present & path_present, c(1, 12, 17)]
  colnames(ko_terms) <- c("Gene", "KO", "Pathway")
  ko_terms$Sample <- sample_id
  
  kegg_list[[sample_id]] <- ko_terms
}

ko_pathway_df <- bind_rows(kegg_list)
write_tsv(ko_pathway_df, "ko_to_pathway_mapping.tsv")
```

Summarise KEGG pathways per sample
```bash
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)

# Assuming your `ko_pathway_df` contains columns: Sample, KEGG_Pathway
# If not loaded yet, run this:
ko_pathway_df <- readr::read_tsv("ko_to_pathway_mapping.tsv")

# Summarize counts per sample and pathway
kegg_pathway_counts <- ko_pathway_df %>%
  group_by(Sample, KEGG_Pathway) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = KEGG_Pathway, values_from = Count, values_fill = 0)

# Save output for reference
readr::write_tsv(kegg_pathway_counts, "kegg_pathway_counts_per_sample.tsv")

# Columns: Gene | KO | Pathway | Sample
ko_pathway_df <- read_tsv("ko_to_pathway_mapping.tsv")

library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(tibble)

# Flatten multiple pathways per KO
kegg_counts <- ko_pathway_df %>%
  separate_rows(Pathway, sep = ",") %>%
  group_by(Sample, Pathway) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Pathway, values_from = Count, values_fill = 0)

# Save for inspection
write.table(kegg_counts, "KEGG_pathway_counts_per_sample.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make sample names rownames
kegg_matrix <- kegg_counts %>%
  column_to_rownames("Sample") %>%
  as.matrix()

# Optional: scale values by row
kegg_scaled <- t(scale(t(kegg_matrix)))

# Plot heatmap
pheatmap(kegg_scaled,
         main = "KEGG Pathway Counts per Sample",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize = 9)

# Top 20 pathways across all samples
top_pathways <- ko_pathway_df %>%
  separate_rows(Pathway, sep = ",") %>%
  count(Pathway, sort = TRUE) %>%
  slice_head(n = 20)

# Filter out broad overview KEGG pathways (optional)
top_pathways_filtered <- top_pathways %>%
  filter(!Pathway %in% c("ko00000", "ko00001", "ko01000"))

# Redraw the plot
ggplot(top_pathways_filtered, aes(x = reorder(Pathway, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top KEGG Pathways (excluding overview categories)", x = "KEGG Pathway", y = "Count") +
  theme_minimal()


