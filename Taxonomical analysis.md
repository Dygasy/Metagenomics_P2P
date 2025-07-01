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

Correlation Analysis (is this the same as maaslin2?)
1. Correlate specific taxa abundances with clinical metadata or other experimental variables

Are we doing this for functional data generated as well?

Use of PICRUSt2 or Tax4Fun (GTDB primarily provides taxonomic classification, but the taxa identified can be mapped indirectly to known genomes, facilitating subsequent functional predictions using tools mentioned)
