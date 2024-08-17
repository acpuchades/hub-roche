library(readr)

library(plotly)
library(htmlwidgets)

variants.pca <- read_tsv("output/common-variants/plink.eigenvec")

plot_ly(
  variants.pca,
  text = ~IID, x = ~PC1, y = ~PC2, z = ~PC3,
  hoverinfo = "text", size = 0.5
) |>
  saveWidget("output/common-variants/variants-pca.html")
