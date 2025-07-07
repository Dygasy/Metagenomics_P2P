MEGAN's gtdbpath_to_percent and gtdbname_to_percent exports has already calculated relative abundance for you. 

MEGAN's gtdbpath_to_percent
Contains the full taxonomy path  (e.g., Bacteria;Firmicutes;Clostridia;...)
Values = percent of reads assigned to that path per sample
Best for Krona or rank-based stacked plots

MEGAN's gtdbname_to_percent
Contains just the taxon name at a specific level (eg: Clostridia, Bacteroides) depending on the classification you're going for
Still shows % abundance per sample
useful for barplots, heatmaps,or clustering at fixed levels

Perform Taxonomical Composition (Stacked bar plots) for all taxa (Phylum --> Species)
R script
```bash
# === Load required libraries ===
library(tidyverse)
library(readr)

# === 1. Load the relative abundance data ===
df <- read_tsv("E:/Krona_results/05_Taxonomy/Relative abundance/Class_relativeabundance.txt")

# === 2. Tidy the data: reshape to long format ===
df_long <- df %>%
  rename(Class = 1) %>%  # Ensure first column is named 'Class'
  pivot_longer(-Class, names_to = "Sample", values_to = "Abundance")

# === 3. Filter to top 20 most abundant classes, group others ===
top_classes <- df_long %>%
  group_by(Class) %>%
  summarise(total = sum(Abundance, na.rm = TRUE)) %>%
  slice_max(total, n = 20) %>%
  pull(Class)

df_long <- df_long %>%
  mutate(Class = ifelse(Class %in% top_classes, Class, "Others"))

# === 4. Plot ===
ggplot(df_long, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal(base_size = 13) +
  labs(title = "Relative Abundance at Class Level",
       y = "Relative abundance (%)", x = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")
```

Perform heatmap using relative abundance table (specific taxa) in R
```bash
# Install once if not yet installed
install.packages("pheatmap")

# Load packages
library(readr)
library(pheatmap)
library(dplyr)

# === 1. Load data ===
df <- read_tsv("E:/Krona_results/05_Taxonomy/Relative abundance/Class_relativeabundance.txt")

# === 2. Convert to matrix ===
df_mat <- df %>%
  rename(Taxon = 1) %>%           # Ensure first column is named "Taxon"
  column_to_rownames("Taxon")     # Use taxa as rownames

# === 3. Convert to numeric matrix ===
df_matrix <- as.matrix(df_mat)

# === 4. Generate heatmap ===
pheatmap(df_matrix,
         scale = "none",                   # no row scaling
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 8,
         fontsize_col = 8,
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Class-level Relative Abundance Heatmap")
```


```bash
# Set factor levels for Sex if needed
alpha_combined$Sex <- factor(alpha_combined$Sex, levels = c("male", "female"))

# Pairwise comparison setup
sex_pairwise <- list(c("male", "female"))

# Updated plotting function to show numeric p-values
plot_metric_by_sex <- function(data, metric, y_label) {
  ggplot(data, aes(x = Sex, y = .data[[metric]], fill = Sex)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Sex"), y = y_label, x = "Sex") +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = sex_pairwise, 
                       label = "p.format")  # this shows actual p-values
}

# Generate plots
plot_sex_chao1   <- plot_metric_by_sex(alpha_combined, "Chao1", "Chao1")
plot_sex_shannon <- plot_metric_by_sex(alpha_combined, "Shannon", "Shannon")
plot_sex_simpson <- plot_metric_by_sex(alpha_combined, "Simpson", "Simpson")
plot_sex_fisher  <- plot_metric_by_sex(alpha_combined, "Fisher", "Fisher")

# Combine and display
ggarrange(plot_sex_chao1, plot_sex_shannon, plot_sex_simpson, plot_sex_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
```
follow suit for the rest of the other metadata factors 


Export absolute and relative abundance GTDB data in different rankings out to generate OTU table in MEGAN community.
R libraries most commonly used: ggplot2, phyloseq, vegan, microboime, complexHeatmap, 
extract from pathcounts - 
Look at Alpha Diversity 
Compute within-sample diversity: consolidate them in a page using R, (vegan, phyloseq2) 
1. Chao1
2. Shannon
3. Fisher
4. Simpson

Look at CLR and rarefaction too. Tutorial on transformations: https://joey711.github.io/phyloseq/
R-Script: 

Look at Beta Diversity
Compute between-sample distances: Bray-Curtis with PERMANOVA, visualise ecological patterns:
1. PCA
2. PCoA
3. NMDS
4. UPGMA tree (Hierarchical clustering- based on pairwise distances can choose Bray-curtis or UniFrac)
5. we did perform adonis - better to expand
5.Perform PERMANOVA - test for significant differences in community composition between groups. 

R-Script:

Look at Differential Abundance Analysis 
Compute and compare taxa between groups (extracted from the different metadata factors available) 
Use: Maaslin2 (PS, no longer have LEfSe)
can consider DESeq2 and ALDEx2 (for compositional data)

R-script:

Look at Taxonomic Barplots and Heatmaps
Aggregate by Species, Genus, Order, Class, Phylum
Taxonomic composition analysis
1. summarise the relative abundance or presence of taxa at different taxonomic levels
2. generate bar plots/stacked bar charts to visualise taxonomic distributions
Heatmaps and Clustering
1. Visualise taxa abundance across samples through heatmaps
2. Perform hierarchical clustering or k-means clustering to identify samples with similar taxonomic profiles

R-script for taxonomical Barplots:
```bash
# === Load libraries ===
library(readxl)
library(tidyverse)
library(pheatmap)

# === Load metadata and sort by BMI ===
metadata <- read_excel("E:/Krona_results/05_Taxonomy/Metadata/Clinical/Meta_Data_Stool_2024-12-06.xlsm")
metadata_sorted <- metadata %>%
  arrange(desc(BMI))  # Replace with correct BMI column name

# === Define taxonomic levels and file paths ===
levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
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
    select(SampleID = `Your_SampleID_Column`, BMI) %>%  # Replace with actual sample ID column
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
```

Perform Cluster Heatmap:
Generating a heatmap of taxa abundance, where samples(rows) are clustered based on a specific variable, and taxa (columns) are the top 20 most abundant (eg; species level). Basically, reordering the rows(samples) in the heatmap based on similarity in taxa abundance, optionally using BMI as an annotation. 

Hierarchical clustering groups samples with similar microbial compositions:
in order to explore natural groupings

R-script: 
```bash
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

# === 2. Load genus-level relative abundance ===
abundance <- read_tsv("E:/Krona_results/05_Taxonomy/Relative abundance/Genus_relativeabundance.txt", show_col_types = FALSE)
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
  main = "Top 20 Genus Abundance (Clustered by Microbiome Profile, Annotated by BMI)",
  scale = "row",
  fontsize_row = 7,
  fontsize_col = 8
)
```

You basically used Euclidean distance between all pairs of samples, based on their scaled (z-score standardised) genus abundances. Apply hierarchical clustering using the complete linkage method. Group similar microbiome profiles together and draw dendogram branches accordingly. 

Perform Species Accumulation Curve: For human gut microbiome, you typically observe 500-2000 species depending on :
1. sequencing depth
2. sample prep
3. Analysis pipeline (taxonomic classifiers & thresholds)
4. Database used (eg; GTDB, RefSeq, MetaPhlAn, etc)
5. For species-level resolution, ~1000 species detected across ~100 samples is realistic and consistent with published studies.

Generate scaftigs (scaffold contigs) is a common step after a metagenomic assembly using tools like MEGAHIT
Use of scaftigs for downstream analyses
filter and analyse only longer assembled sequences
feed into gene prediction, binning, or functional annotation tools


```bash
# === 1. Load required libraries ===
library(tidyverse)
library(vegan)

# === 2. Load absolute species abundance table ===
# Replace with your actual file path
df <- read_tsv("E:/Krona_results/05_Taxonomy/Absolute Abundance/Species_abundance.txt")

# === 3. Convert to matrix: rows = species, columns = samples ===
df_matrix <- df %>%
  column_to_rownames(var = colnames(df)[1]) %>%
  as.matrix()

# === 4. Transpose: now rows = samples, columns = species ===
df_matrix <- t(df_matrix)

# === 5. Compute species accumulation curve ===
spec_acc <- specaccum(df_matrix, method = "random")

# === 6. Plot the accumulation curve ===
plot(
  spec_acc,
  xlab = "Number of samples",
  ylab = "Accumulated observed species",
  main = "Species Accumulation Curve",
  col = "blue",
  ci.type = "polygon",  # confidence interval shading
  ci.col = "lightblue",
  ci.lty = 0,
  lwd = 2
)
```
X-axis (Number of samples): As you add more samples, you're introducing more diversity and potentially detecting new species
Y-axis (Accumulated observed species): Total number of unique species detected across all samples up to that point. 
Blue curve: Average number of observed species per number per samples
shaded area: Confidence interval from random permutations (species accumulation depends on sample order)


R-script: 
```bash
# === Load required libraries ===
library(readxl)
library(tidyverse)
library(pheatmap)

# === Load metadata and abundance ===
metadata <- read_excel("E:/Krona_results/05_Taxonomy/Metadata/Clinical/Meta_Data_Stool_2024-12-06.xlsm")
abundance <- read_tsv("E:/Krona_results/05_Taxonomy/Relative abundance/Genus_relativeabundance.txt", show_col_types = FALSE)

# === Rename first column ===
colnames(abundance)[1] <- "Taxon"

# === Select top 20 most abundant genera ===
top_genera <- abundance %>%
  pivot_longer(-Taxon, names_to = "Sample", values_to = "Abundance") %>%
  group_by(Taxon) %>%
  summarise(total = sum(Abundance, na.rm = TRUE)) %>%
  slice_max(total, n = 20) %>%
  pull(Taxon)

abundance_top <- abundance %>%
  filter(Taxon %in% top_genera)

# === Transpose abundance to wide format ===
abund_matrix <- abundance_top %>%
  pivot_longer(-Taxon, names_to = "SampleID", values_to = "Abundance") %>%
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0)

# === Merge with BMI metadata ===
merged <- metadata %>%
  rename(SampleID = `Your_SampleID_Column`) %>%  # ⛔ Replace this
  select(SampleID, BMI) %>%
  inner_join(abund_matrix, by = "SampleID") %>%
  drop_na(BMI)

# === Prepare matrix for heatmap ===
rownames_mat <- merged$SampleID
annotation_row <- merged %>% select(BMI)
heatmap_data <- merged %>% select(-SampleID, -BMI)
rownames(heatmap_data) <- rownames_mat

# === Generate heatmap clustered by BMI ===
pheatmap(
  mat = as.matrix(heatmap_data),
  annotation_row = annotation_row,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  main = "Top 20 Genus Abundance (Clustered by BMI)",
  scale = "row",
  fontsize_row = 7,
  fontsize_col = 8
)
```


Correlation Analysis (is this the same as maaslin2?)
1. Correlate specific taxa abundances with clinical metadata or other experimental variables

Use actual measurements (Continuous Values)
Good when you want to:
1. Explore correlations (eg; Spearman/Pearson correlation between microbial abundance and BMI/LDL)
2. Visualise trends across gradients (eg; Bacteroides increases with higher BMI)
3. When your sample size is small (categorising it might oversimplify)
4. When planning to do correlation heatmaps or regression
5. for heatmaps, start with actual values, it's more data-rich

Use categories (Desirable/High/Obese) 
Good when you want to:
1. Perform comparative statistics (eg: Kruskal-Wallis or ANOVA between groups)
2. Cluster samples by condition or group in heatmaps
3. Simplify interpretation for presentation or publication
4. When you want clear visual stratification
5. aim to highlight group-level differences (Obese vs lean microbiota profiles)
6. Suitable for summary barplots, statistical comparisons, or LEfSe-like analysis

Are we doing this for functional data generated as well?

Use of PICRUSt2 or Tax4Fun (GTDB primarily provides taxonomic classification, but the taxa identified can be mapped indirectly to known genomes, facilitating subsequent functional predictions using tools mentioned)
