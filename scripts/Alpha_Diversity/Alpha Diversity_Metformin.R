# Categorize by Metformin use
alpha_combined$Metformin_cat <- factor(
  alpha_combined$Metformin,
  levels = c("no", "yes")
)

# Pairwise comparisons (only 2 groups here: yes vs no)
metformin_pairs <- list(c("no", "yes"))

# Plotting function
plot_metric_by_metformin <- function(data, metric, y_label) {
  ggplot(data, aes(x = Metformin_cat, y = .data[[metric]], fill = Metformin_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Metformin Use"), x = "Metformin", y = y_label) +
    stat_compare_means(method = "wilcox.test") +
    stat_compare_means(comparisons = metformin_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_metformin_chao1   <- plot_metric_by_metformin(alpha_combined, "Chao1", "Chao1")
plot_metformin_shannon <- plot_metric_by_metformin(alpha_combined, "Shannon", "Shannon")
plot_metformin_simpson <- plot_metric_by_metformin(alpha_combined, "Simpson", "Simpson")
plot_metformin_fisher  <- plot_metric_by_metformin(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_metformin_chao1, plot_metformin_shannon, plot_metformin_simpson, plot_metformin_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
