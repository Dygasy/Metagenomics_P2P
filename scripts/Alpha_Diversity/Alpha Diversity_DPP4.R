# Categorize by DPP4 inhibitor use
alpha_combined$DPP4_cat <- factor(
  alpha_combined$`DPP4 inhibitor`,
  levels = c("no", "yes")
)

# Pairwise comparisons
dpp4_pairs <- list(c("no", "yes"))

# Plotting function
plot_metric_by_dpp4 <- function(data, metric, y_label) {
  ggplot(data, aes(x = DPP4_cat, y = .data[[metric]], fill = DPP4_cat)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by DPP4 Inhibitor Use"), x = "DPP4 Inhibitor", y = y_label) +
    stat_compare_means(method = "wilcox.test") +
    stat_compare_means(comparisons = dpp4_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate plots
plot_dpp4_chao1   <- plot_metric_by_dpp4(alpha_combined, "Chao1", "Chao1")
plot_dpp4_shannon <- plot_metric_by_dpp4(alpha_combined, "Shannon", "Shannon")
plot_dpp4_simpson <- plot_metric_by_dpp4(alpha_combined, "Simpson", "Simpson")
plot_dpp4_fisher  <- plot_metric_by_dpp4(alpha_combined, "Fisher", "Fisher")

# Display all plots
ggarrange(plot_dpp4_chao1, plot_dpp4_shannon, plot_dpp4_simpson, plot_dpp4_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
