library(tidyverse)
library(gridExtra)
library(brms)
library(bayestestR)
library(tidybayes)
library(ggpubr)

# Read in data file
mod_dataR <- read_csv("Model_Data.csv")

# Model for nest site selection
mod1 <- brm(Nest_Random ~ tph*Died_Nesting + gc*Died_Nesting + cc*Died_Nesting +
            maxveg*Died_Nesting + voveg*Died_Nesting + (1|Site),
            data = mod_dataR, chains = 4, cores = 4, family = bernoulli(), iter = 5000,
            warmup = 2000)

print(summary(mod1), digits = 3)

# Extract plots for gc, voveg, and tph
effects <- conditional_effects(mod1, effects = c("gc", "voveg", "tph"))

# Plot function
plot_effect <- function(effect_data, h, x_var, x_label, y_label = NULL, hist_height = 0.05) {
  h <- h %>%
    mutate(pct_adj = ifelse(Nest_Random == 0, pct * hist_height, 1 - pct * hist_height))
  
  ggplot() +
    geom_segment(data = h, size = 3, show.legend = FALSE,  # Shrink the size of histograms
                 aes_string(x = x_var, xend = x_var, y = "Nest_Random", yend = "pct_adj", colour = "factor(Nest_Random)")) +
    geom_ribbon(data = effect_data, aes_string(x = x_var, ymin = "lower__", ymax = "upper__"), fill = "grey80", alpha = 0.5) +
    geom_line(data = effect_data, aes_string(x = x_var, y = "estimate__"), color = "grey50", size = 1) +
    scale_y_continuous(limits = c(-0.02, 1.02)) +
    scale_x_continuous(name = x_label) +
    labs(y = y_label) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
}

# Summarize data for tph, gc, and voveg
h_tph <- summarize_data(mod_dataR_clean, "tph")
h_gc <- summarize_data(mod_dataR_clean, "gc")
h_voveg <- summarize_data(mod_dataR_clean, "voveg")

# Create plots
plot_tph <- plot_effect(effects[["tph"]], h_tph, "tph", "Trees per ha", hist_height = 0.85)
plot_gc <- plot_effect(effects[["gc"]], h_gc, "gc", "Ground cover (%)", "Probability of selection", hist_height = 0.85)
plot_voveg <- plot_effect(effects[["voveg"]], h_voveg, "voveg", "Visual obstruction (cm)", hist_height = 0.85)

# Combine plots
combined_plot <- (plot_gc | plot_voveg | plot_tph) +
  plot_layout(ncol = 3, widths = c(1, 1, 1))

# Print the final combined plot
print(combined_plot)

# Save the combined plot at 600 dpi
ggsave("Figure1_Updated.png", plot = combined_plot, dpi = 600, width = 12, height = 4, units = "in")
