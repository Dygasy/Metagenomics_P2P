flowchart LR
    T1["Input: GTDB counts"]
    T2["Tidy & filter data"]
    T3["ggplot2 stacked barplots"]
    T4["pheatmap heatmaps"]

    A1["Compute diversity indices"]
    A2["Wilcoxon tests"]
    A3["PERMANOVA"]
    A4["NMDS / PCoA"]

    D1["Prepare OTU & metadata"]
    D2["Maaslin2"]
    D3["Significant taxa"]

    S1["Species accumulation curve"]

    %% Connections
    T1 --> T2 --> T3
    T2 --> T4
    T2 --> A1 --> A2
    A1 --> A3 --> A4
    T2 --> D1 --> D2 --> D3
    T2 --> S1

    %% Style boxes lightly
    classDef output fill:#fdf6e3,stroke:#000,stroke-width:1px,font-weight:bold;
    class T3,T4,A2,A4,D3,S1 output;

