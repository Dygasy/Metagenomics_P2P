library(ggplot2)
library(ggpubr)

# Set factor levels for Sex if needed
alpha_combined$Sex <- factor(alpha_combined$Sex, levels = c("male", "female"))

# Optional: Pairwise comparison setup
sex_pairwise <- list(c("male", "female"))

# Updated plotting function with violin + boxplot + jitter
plot_metric_by_sex <- function(data, metric, y_label) {
  ggplot(data, aes(x = Sex, y = .data[[metric]], fill = Sex)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    theme_minimal() +
    labs(title = paste(y_label, "Diversity by Sex"), y = y_label, x = "Sex") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = sex_pairwise)
}

# Now generate the plots
plot_sex_chao1   <- plot_metric_by_sex(alpha_combined, "Chao1", "Chao1")
plot_sex_shannon <- plot_metric_by_sex(alpha_combined, "Shannon", "Shannon")
plot_sex_simpson <- plot_metric_by_sex(alpha_combined, "Simpson", "Simpson")
plot_sex_fisher  <- plot_metric_by_sex(alpha_combined, "Fisher", "Fisher")

# Combine and display
ggarrange(plot_sex_chao1, plot_sex_shannon, plot_sex_simpson, plot_sex_fisher,
          ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
