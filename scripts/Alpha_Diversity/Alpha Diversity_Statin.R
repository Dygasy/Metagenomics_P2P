# Categorize Statin use (already yes/no, but we format as factor)
alpha_combined$Statin_cat <- factor(
  alpha_combined$Statin,
  levels = c("no", "yes"),
  labels = c("No Statin", "Statin")
)

# Pairwise comparison group
statin_pairs <- combn(levels(alpha_combined$Statin_cat), 2, simplify = FALSE)

# Plotting function for Statin
plot_metric_by_statin <- function(data, metric, y_label) {
  ggplot(data, aes(x = Statin_cat, y = .data[[metric]], fill = Statin_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Statin Use"), x = "Statin Use", y = y_label) +
    stat_compare_means(method = "wilcox.test") +
    stat_compare_means(comparisons = statin_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_statin_chao1   <- plot_metric_by_statin(alpha_combined, "Chao1", "Chao1")
plot_statin_shannon <- plot_metric_by_statin(alpha_combined, "Shannon", "Shannon")
plot_statin_simpson <- plot_metric_by_statin(alpha_combined, "Simpson", "Simpson")
plot_statin_fisher  <- plot_metric_by_statin(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_statin_chao1, plot_statin_shannon, plot_statin_simpson, plot_statin_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
