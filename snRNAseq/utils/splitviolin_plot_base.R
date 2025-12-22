library(tools)
library(Seurat)
library(ggplot2)
library(introdataviz)
library(forcats)
library(ggpubr)

# Function to add module scores
add_module_score <- function(seurat_obj, gene_list_path, module_name) {
  gene_list <- read.csv(gene_list_path)
  gene_list <- toTitleCase(tolower(gene_list$Gene.ID))
  features <- rownames(seurat_obj@assays$originalexp)
  
  matching_genes <- gene_list %in% features
  cat("Number of matching genes for", module_name, ":", sum(matching_genes), "\n")
  
  gene_list_filtered <- gene_list[matching_genes]
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = list(gene_list_filtered),
    ctrl = 500,
    name = module_name
  )
  return(seurat_obj)
}

my_palette <- function(n) {
  colors <- c("#7A7A7A", "#DC9F41", "#4DAF4A", "#984EA3", "#FF7F00")  # Custom colors
  return(colors[1:n])
}
## Make split violin plot
plot_split_violin <- function(plot_data, tglname, x_var, y_var, tsuffix = 'test',fsizew=8,fsizeh=4) {
  p <- ggplot(plot_data, aes(
    x = .data[[x_var]],
    y = .data[[tglname]],
    fill = .data[[y_var]]
  )) +
    geom_split_violin(width = 1) +
    scale_fill_manual(values = my_palette(length(unique(plot_data[[y_var]]))))+
    geom_point(aes(group = .data[[y_var]]),
               position = position_jitterdodge(jitter.width = 0.4),
               size = 0.01,
               alpha = 0.5) +
    theme(
      panel.grid.major = element_blank(),  # Removes major grid lines
      panel.grid.minor = element_blank(),  # Removes minor grid lines
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),  # Adds black axis lines
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(y = tglname)  # Dynamically set y-axis text to the gene name

  # Add p-value using ggpubr
  p <- p + stat_compare_means(
    aes(label = after_stat(p.format)),  # Display p-value
    method = "wilcox.test",             # Statistical test (e.g., Wilcoxon rank-sum test)
    label.x = 1.5,                      # Position of the p-value on the x-axis
    label.y = max(plot_data[[tglname]], na.rm = TRUE) * 1.1  # Position on the y-axis
  )

  # Save the plot
  ggsave(
    filename = paste0("figures/TSplit_ViolinPlot_", tglname, tsuffix, ".pdf"),
    plot = p,
    width = fsizew,
    height = fsizeh
  )
}
## make statitis about numbers and expression level
summarize_group_stats <- function(plot_data, tglname, x_var, y_var, tsuffix = 'test') {
  require(dplyr)
  
  # Create summary statistics
  summary_df <- plot_data %>%
    group_by(.data[[x_var]], .data[[y_var]]) %>%
    summarise(
      n_cells = sum(!is.na(.data[[tglname]])),
      mean_expression = mean(.data[[tglname]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(
      time_point = 1,
      condition = 2
    )
  
  # Create filename with tsuffix
  filename <- paste0("csv/",tglname, "_group_stats_", tsuffix, ".csv")
  
  # Save to CSV
  write.csv(summary_df, file = filename, row.names = FALSE)
  
  message("Saved group statistics to: ", filename)
  
  return(summary_df)
}

summarize_group_stats <- function(plot_data, tglname, x_var, y_var, tsuffix = 'test') {
  require(dplyr)
  
  # Calculate summary statistics
  summary_df <- plot_data %>%
    group_by(across(c(.data[[x_var]], .data[[y_var]]))) %>%
    summarise(
      n_cells = sum(!is.na(.data[[tglname]])),
      mean_expression = mean(.data[[tglname]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(time_point = 1, condition = 2)
  
  # Calculate Wilcoxon p-values for each time point
  pval_df <- plot_data %>%
    group_by(across(.data[[x_var]])) %>%
    summarise(
      p_value = tryCatch({
        wilcox.test(.data[[tglname]] ~ .data[[y_var]], exact = FALSE)$p.value
      }, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    rename(time_point = 1)
  
  # Merge statistics with p-values
  summary_df <- summary_df %>%
    left_join(pval_df, by = "time_point") %>%
    mutate(p_value = formatC(p_value, format = "e", digits = 3))
  
  # Save to CSV
  filename <- paste0("csv/",tglname, "_group_stats_", tsuffix, ".csv")
  write.csv(summary_df, file = filename, row.names = FALSE)
  message("Saved group statistics to: ", filename)
  
  return(summary_df)
}