# Rscript docs/benchmarks/gpu-sketch-geuvadis/plot.R  (from the repository root)
# Rebuilds gpu-sketch-geuvadis.png from gpu-sketch-geuvadis.csv.
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })
dir <- "docs/benchmarks/gpu-sketch-geuvadis"
d <- read.csv(file.path(dir, "gpu-sketch-geuvadis.csv"))
cols <- c("GPU lookup (Metal)" = "#D55E00", "CPU, 16 threads (--no-gpu)" = "#0072B2")
curve <- subset(d, sample != "ERR188428")
best <- do.call(rbind, lapply(split(curve, list(curve$pairs, curve$v), drop = TRUE),
                              function(g) g[which.min(g$wall), ]))
panel <- function(y, title) {
  ggplot(curve, aes(pairs, .data[[y]], colour = v)) +
    geom_point(alpha = 0.6, size = 1.8) +
    geom_line(data = best[order(best$pairs), ], linewidth = 0.8) +
    scale_colour_manual(values = cols, name = NULL) +
    expand_limits(x = 0, y = 0) +
    labs(x = "Read pairs (millions), ERR188044 and its subsamples", y = "Seconds", title = title) +
    theme_cowplot(11) + theme(legend.position = "none")
}
full <- subset(d, sample != "ERR188044 subsample")
full_best <- do.call(rbind, lapply(split(full, list(full$sample, full$v), drop = TRUE),
                                   function(g) g[which.min(g$wall), ]))
c_panel <- ggplot(full_best, aes(sample, wall, fill = v)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f s", wall)), position = position_dodge(0.8),
            vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = cols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Seconds", title = "Whole run, full samples") +
  theme_cowplot(11) + theme(legend.position = "none")
legend <- get_plot_component(panel("wall", "") + theme(legend.position = "bottom"),
                             "guide-box-bottom", return_all = TRUE)
grid <- plot_grid(panel("wall", "Whole run"), panel("mapping", "Mapping phase"), c_panel,
                  nrow = 1, labels = c("A", "B", "C"), rel_widths = c(1, 1, 0.8))
ggsave(file.path(dir, "gpu-sketch-geuvadis.png"),
       plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.08)),
       width = 11, height = 4.2, dpi = 200, bg = "white")
