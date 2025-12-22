library(ggplot2)
library(reshape2)

generate_boxplot <- function(tdf, output_pdf,fig_width=8, fig_height=6) {
  df_melted <- melt(tdf, id.vars = 'X')
  X <- as.matrix(tdf[2:nrow(tdf), 2:ncol(tdf)])
  tname <- names(sort(colMeans(X, na.rm = TRUE), decreasing = TRUE))
  df_melted$variable <- factor(df_melted$variable, levels = as.list(tname))
  pdf(output_pdf, width = fig_width, height = fig_height)
  p <- ggplot(df_melted, aes(x = factor(variable), y = value, fill = factor(variable))) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 16,
                 outlier.size = 2, notch = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  print(p)
  dev.off()
}