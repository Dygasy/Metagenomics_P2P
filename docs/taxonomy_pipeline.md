# Taxonomy Analysis Pipeline
```mermaid
flowchart LR
 subgraph Taxonomic Composition ["Taxonomic Composition Analysis"]
    T1["Input: GTDB gtdbpath_to_count &amp; gtdbname_to_count"]
    T2["Tidy & filter data (top 20 taxa)"]
    T3["ggplot2 stacked barplots"]
    T4["pheatmap cluster heatmaps"]

 subgraph Alpha_Beta_Diversity["Alpha & Beta Diversity Analysis"]
    A1["Compute diversity indices (Chao1, Shannon, Simpson, Fisher) via vegan"]
    A2["Wilcoxon tests (sex, BMI, etc"]
    A3["PERMANOVA (adonis2) on Bray Curtis"]
    A4["NMDS / PCoA ordinations"]

 subgraph Differential_Abundance["Differential Abundance"]
    D1["Prepare OTU & metadata tables"]
    D2["Maaslin2 linear models"]
    D3["Output: Significant taxa (FDR&lt;0.05)"]

 subgraph Species_Accumulation["Species Accumulation"]
    S1["vegan specaccum curve"]

    %% Connections
    T1 --> T2 & A1 & D1 & S1
    T2 --> T3 & T4
    A1 --> A2 & A3
    A3 --> A4
    D1 --> D2
    D2 --> D3

    %% Style boxes lightly
    classDef output fill:#fff2cc,stroke:#000,stroke-width:2px,font-weight:bold;
    class T3,T4,A2,A4,D3,S1 output;

