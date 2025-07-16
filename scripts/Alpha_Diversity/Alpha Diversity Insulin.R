# Categorize by Insulin use
alpha_combined$Insulin_cat <- factor(
  alpha_combined$Insulin,
  levels = c("no", "yes")
)

# Pairwise comparisons (2 groups)
insulin_pairs <- list(c("no", "yes"))

# Plotting function
plot_metric_by_insulin <- function(data, metric, y_label) {
  ggplot(data, aes(x = Insulin_cat, y = .data[[metric]], fill = Insulin_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Insulin Use"), x = "Insulin", y = y_label) +
    stat_compare_means(method = "wilcox.test") +
    stat_compare_means(comparisons = insulin_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_insulin_chao1   <- plot_metric_by_insulin(alpha_combined, "Chao1", "Chao1")
plot_insulin_shannon <- plot_metric_by_insulin(alpha_combined, "Shannon", "Shannon")
plot_insulin_simpson <- plot_metric_by_insulin(alpha_combined, "Simpson", "Simpson")
plot_insulin_fisher  <- plot_metric_by_insulin(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_insulin_chao1, plot_insulin_shannon, plot_insulin_simpson, plot_insulin_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
