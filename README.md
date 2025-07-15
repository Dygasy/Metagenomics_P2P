# Metagenomics_V2

This repository contains a comprehensive pipeline for metagenomic analyses developed for our collaborative study, encompassing both taxonomic and functional aspects:


Features

Taxonomic Composition
- **Visualisation:** stacked barplots, heatmaps (ggplot2, pheatmap)
- **Alpha Diversity:** Chao1, Shannon, Simpson, Fisher indices (vegan)
- **Beta Diversity:** Bray-Curtis dissimilarityy, NMDS ordinations, PERMANOVA
- **Differential Abundance:** Maaslin2 linear models to identify taxa by metadata 

Functional Profiling

- **Functional Potential:** Annotates and quantifies metabolic pathways and gene families (using tools like eggNOG-mapper)
- **Visualisation:** Pathway coverage and abundance barplots, KEGG pathway maps
- **Comparative Functional Analysis:** PERMANOVA or differential pathway abundance testing. 



---

## Repository structure


```
data/     # input datasets (GTDB outputs, eggNOG-mapper outputs,metadata files)
docs/     # detailed markdown workflow notes, Mermaid crafted diagrams
figures/  # generated plots (barplots, heatmaps, NMDS, functional profiles)
scripst/  # bash, Python and R scripts for processing, analysis and visualisation
```


--- 


## Quickstart

Clone the repository:


```bash
git clone https://github.com/Dygasy/Metagenomics_V2.git
cd Metagenomics_V2
```


Create your environment (example if using miniconda):

```bash
conda env create -f environment.yml
conda activate metagenomics
```

Then run analyses, eg.:

```bash
Rscript scripts/taxonomic_analysis.R
```

--- 

## license
This project is licensed under the MIT license


---

## Citation

If you use this pipeline or adapt it in your research, please cite this repository (or the forthcoming manuscript).

---

## Contact

For questions, open an issue or contact [Dygasy on Github](https://github.com/Dygasy).
