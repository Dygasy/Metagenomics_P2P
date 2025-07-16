# Categorize Age directly with specified labels
alpha_combined$Age_group <- cut(
  alpha_combined$Age,
  breaks = c(17, 39, 59, Inf),
  labels = c("18–39", "40–59", "60 and above"),
  right = TRUE
)

# Set order of age groups
alpha_combined$Age_group <- factor(alpha_combined$Age_group,
                                   levels = c("18–39", "40–59", "60 and above"))

# Define pairwise comparisons
age_pairs <- combn(levels(alpha_combined$Age_group), 2, simplify = FALSE)

# Plotting function
plot_metric_by_age <- function(data, metric, y_label) {
  ggplot(data, aes(x = Age_group, y = .data[[metric]], fill = Age_group)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Age Group"), x = "Age Group", y = y_label) +
    stat_compare_means(method = "kruskal.test") +
    stat_compare_means(comparisons = age_pairs, method = "wilcox.test", label = "p.signif")
}

# Generate and display plots
plot_age_chao1   <- plot_metric_by_age(alpha_combined, "Chao1", "Chao1")
plot_age_shannon <- plot_metric_by_age(alpha_combined, "Shannon", "Shannon")
plot_age_simpson <- plot_metric_by_age(alpha_combined, "Simpson", "Simpson")
plot_age_fisher  <- plot_metric_by_age(alpha_combined, "Fisher", "Fisher")

ggarrange(plot_age_chao1, plot_age_shannon, plot_age_simpson, plot_age_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
