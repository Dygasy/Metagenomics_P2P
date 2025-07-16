# Set factor levels (optional, to control plot order)
alpha_combined$Ethnicity <- factor(alpha_combined$Ethnicity,
                                   levels = c("chinese", "malay", "indian", "others"))

# Generate all pairwise comparisons for Ethnicity groups
ethnicity_levels <- levels(alpha_combined$Ethnicity)
ethnicity_pairs <- combn(ethnicity_levels, 2, simplify = FALSE)

# Plotting function
plot_metric_by_ethnicity <- function(data, metric, y_label) {
  ggplot(data, aes(x = Ethnicity, y = .data[[metric]], fill = Ethnicity)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Ethnicity"), y = y_label, x = "Ethnicity") +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = ethnicity_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_eth_chao1   <- plot_metric_by_ethnicity(alpha_combined, "Chao1", "Chao1")
plot_eth_shannon <- plot_metric_by_ethnicity(alpha_combined, "Shannon", "Shannon")
plot_eth_simpson <- plot_metric_by_ethnicity(alpha_combined, "Simpson", "Simpson")
plot_eth_fisher  <- plot_metric_by_ethnicity(alpha_combined, "Fisher", "Fisher")

# Combine and show plots
ggarrange(plot_eth_chao1, plot_eth_shannon, plot_eth_simpson, plot_eth_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
