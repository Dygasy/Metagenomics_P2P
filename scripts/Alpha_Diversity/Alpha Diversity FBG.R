# Categorize FBG
alpha_combined$FBG_cat <- cut(
  alpha_combined$FBG,
  breaks = c(-Inf, 5.9, 6.9, Inf),
  labels = c("Normal", "Pre-diabetes", "Diabetes"),
  right = TRUE
)

# Reorder factor levels
alpha_combined$FBG_cat <- factor(
  alpha_combined$FBG_cat,
  levels = c("Normal", "Pre-diabetes", "Diabetes")
)

# Pairwise comparison groups
fbg_pairs <- combn(levels(alpha_combined$FBG_cat), 2, simplify = FALSE)

# Plotting function
plot_metric_by_fbg <- function(data, metric, y_label) {
  ggplot(data, aes(x = FBG_cat, y = .data[[metric]], fill = FBG_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by FBG"), x = "FBG Category", y = y_label) +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = fbg_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_fbg_chao1   <- plot_metric_by_fbg(alpha_combined, "Chao1", "Chao1")
plot_fbg_shannon <- plot_metric_by_fbg(alpha_combined, "Shannon", "Shannon")
plot_fbg_simpson <- plot_metric_by_fbg(alpha_combined, "Simpson", "Simpson")
plot_fbg_fisher  <- plot_metric_by_fbg(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_fbg_chao1, plot_fbg_shannon, plot_fbg_simpson, plot_fbg_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
