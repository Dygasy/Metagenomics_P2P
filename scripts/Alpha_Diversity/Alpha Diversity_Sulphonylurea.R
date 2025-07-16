# Categorize by Sulphonylurea use
alpha_combined$Sulphonylurea_cat <- factor(
  alpha_combined$Sulphonylurea,
  levels = c("no", "yes")
)

# Pairwise comparisons (2 groups)
sulphonylurea_pairs <- list(c("no", "yes"))

# Plotting function
plot_metric_by_sulphonylurea <- function(data, metric, y_label) {
  ggplot(data, aes(x = Sulphonylurea_cat, y = .data[[metric]], fill = Sulphonylurea_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Sulphonylurea Use"), x = "Sulphonylurea", y = y_label) +
    stat_compare_means(method = "wilcox.test") +
    stat_compare_means(comparisons = sulphonylurea_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_sulphonylurea_chao1   <- plot_metric_by_sulphonylurea(alpha_combined, "Chao1", "Chao1")
plot_sulphonylurea_shannon <- plot_metric_by_sulphonylurea(alpha_combined, "Shannon", "Shannon")
plot_sulphonylurea_simpson <- plot_metric_by_sulphonylurea(alpha_combined, "Simpson", "Simpson")
plot_sulphonylurea_fisher  <- plot_metric_by_sulphonylurea(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_sulphonylurea_chao1, plot_sulphonylurea_shannon,
          plot_sulphonylurea_simpson, plot_sulphonylurea_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
