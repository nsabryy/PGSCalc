library(ggplot2)
library(ggrepel) 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript ks_plot.r <input_file> <output_file>")
}

input_fp <- args[1]
output_fp <- args[2]
fp_label <- sub("_ks\\.csv$", "", basename(input_fp))

ks_results <- read.delim(input_fp)

plot <- ggplot(ks_results, aes(x = `KS.Statistic`, y = `p.value`, label = `Population.Pair`)) +
  geom_point() +
  geom_text_repel() +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "KS Test Statistics vs p-value",
       subtitle = fp_label,
       x = "KS Statistic",
       y = "p-value")

ggsave(filename = output_fp, plot = plot, width = 10, height = 6, dpi = 300)
