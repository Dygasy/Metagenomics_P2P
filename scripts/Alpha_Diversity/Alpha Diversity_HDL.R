# Categorize HDL levels
alpha_combined$HDL_cat <- cut(
  alpha_combined$HDL,
  breaks = c(-Inf, 0.99, 1.5, Inf),
  labels = c("Low", "Desirable", "Optimal"),
  right = TRUE
)

# Set factor level order
alpha_combined$HDL_cat <- factor(alpha_combined$HDL_cat,
                                 levels = c("Low", "Desirable", "Optimal"))

# Pairwise comparisons for HDL
hdl_levels <- levels(alpha_combined$HDL_cat)
hdl_pairs <- combn(hdl_levels, 2, simplify = FALSE)

# Plotting function
plot_metric_by_hdl <- function(data, metric, y_label) {
  ggplot(data, aes(x = HDL_cat, y = .data[[metric]], fill = HDL_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by HDL Category"), x = "HDL Category", y = y_label) +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = hdl_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_hdl_chao1   <- plot_metric_by_hdl(alpha_combined, "Chao1", "Chao1")
plot_hdl_shannon <- plot_metric_by_hdl(alpha_combined, "Shannon", "Shannon")
plot_hdl_simpson <- plot_metric_by_hdl(alpha_combined, "Simpson", "Simpson")
plot_hdl_fisher  <- plot_metric_by_hdl(alpha_combined, "Fisher", "Fisher")

# Combine plots
ggarrange(plot_hdl_chao1, plot_hdl_shannon, plot_hdl_simpson, plot_hdl_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
