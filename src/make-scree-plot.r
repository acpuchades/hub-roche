library(dplyr)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
input_path <- coalesce(args[1], "output/variant-analysis/samples.eigenval")
output_path <- coalesce(args[2], "output/variant-analysis/samples.eigenval.scree.png")

eigenvalues <- read.table(input_path, header = FALSE, col.names = "eigenval") |>
  mutate(
    pc = as.factor(row_number()),
    explained_var = eigenval / sum(eigenval)
  )

ggplot(eigenvalues, aes(pc, explained_var, fill = cumsum(explained_var) < 0.9)) +
  geom_col() +
  geom_point() +
  geom_line(group = 1) +
  labs(x = "PC", y = "% Explained Variance") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  theme(legend.position = "none")
ggsave(output_path, width = 6, height = 4, dpi = 300)
