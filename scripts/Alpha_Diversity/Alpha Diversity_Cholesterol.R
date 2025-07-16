# Categorize Total Cholesterol
alpha_combined$Cholesterol_cat <- cut(
  alpha_combined$Cholesterol,
  breaks = c(-Inf, 5.1, 6.1, Inf),
  labels = c("Desirable", "Borderline High", "High"),
  right = TRUE
)

# Reorder factor levels
alpha_combined$Cholesterol_cat <- factor(
  alpha_combined$Cholesterol_cat,
  levels = c("Desirable", "Borderline High", "High")
)

# Pairwise comparison groups
chol_pairs <- combn(levels(alpha_combined$Cholesterol_cat), 2, simplify = FALSE)

# Plot function
plot_metric_by_chol <- function(data, metric, y_label) {
  ggplot(data, aes(x = Cholesterol_cat, y = .data[[metric]], fill = Cholesterol_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Cholesterol"), x = "Cholesterol Category", y = y_label) +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = chol_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_chol_chao1   <- plot_metric_by_chol(alpha_combined, "Chao1", "Chao1")
plot_chol_shannon <- plot_metric_by_chol(alpha_combined, "Shannon", "Shannon")
plot_chol_simpson <- plot_metric_by_chol(alpha_combined, "Simpson", "Simpson")
plot_chol_fisher  <- plot_metric_by_chol(alpha_combined, "Fisher", "Fisher")

# Combine and display
ggarrange(plot_chol_chao1, plot_chol_shannon, plot_chol_simpson, plot_chol_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
