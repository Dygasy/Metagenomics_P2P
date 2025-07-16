# Categorize HbA1c based on your criteria
alpha_combined$HbA1c_cat <- cut(
  alpha_combined$HbA1c,
  breaks = c(-Inf, 5.9, 6.4, Inf),
  labels = c("Normal", "Prediabetes", "Diabetes"),
  right = TRUE
)

# Set factor order
alpha_combined$HbA1c_cat <- factor(alpha_combined$HbA1c_cat,
                                   levels = c("Normal", "Prediabetes", "Diabetes"))

# Create pairwise combinations for comparisons
hba1c_pairs <- combn(levels(alpha_combined$HbA1c_cat), 2, simplify = FALSE)

# Plotting function by HbA1c category
plot_metric_by_hba1c <- function(data, metric, y_label) {
  ggplot(data, aes(x = HbA1c_cat, y = .data[[metric]], fill = HbA1c_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by HbA1c Category"),
         x = "HbA1c Category", y = y_label) +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = hba1c_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_hba1c_chao1   <- plot_metric_by_hba1c(alpha_combined, "Chao1", "Chao1")
plot_hba1c_shannon <- plot_metric_by_hba1c(alpha_combined, "Shannon", "Shannon")
plot_hba1c_simpson <- plot_metric_by_hba1c(alpha_combined, "Simpson", "Simpson")
plot_hba1c_fisher  <- plot_metric_by_hba1c(alpha_combined, "Fisher", "Fisher")

# Combine all
ggarrange(plot_hba1c_chao1, plot_hba1c_shannon, plot_hba1c_simpson, plot_hba1c_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
