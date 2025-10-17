benchmark_mse <- function(deconvolution, groundtruth, cells_extra = NULL,
                          scatter = TRUE, plot = FALSE, file_name = NULL,
                          width = 16, height = 8) {
  groundtruth = groundtruth[rownames(deconvolution), ] #ensure to have the same order of samples
  
  cell_types = c("B.cells", "B.naive.cells", "B.memory.cells",
                 "Macrophages.cells", "Macrophages.M0", "Macrophages.M1",
                 "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells",
                 "NK.activated", "NK.resting", "NKT.cells", "CD4.cells",
                 "CD4.memory.activated", "CD4.memory.resting", "CD4.naive",
                 "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory",
                 "T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells",
                 "Dendritic.activated.cells", "Dendritic.resting.cells",
                 "Cancer", "Endothelial", "Eosinophils", "Plasma", "Myocytes",
                 "Fibroblasts", "Mast.cells", "Mast.activated.cells",
                 "Mast.resting.cells", "CAF")
  cell_types = c(cell_types, cells_extra)
  
  pattern <- paste0("(_", gsub("\\.", "\\\\.", cell_types), ")$", collapse = "|")
  deconvolution_combinations <- unique(gsub(pattern, "", colnames(deconvolution)))
  deconvolution_combinations = gsub("(BPRNACan3DProMet|BPRNACanProMet|BPRNACan)",
                                    "\\1_", deconvolution_combinations)
  
  scatter_plots = function(deconv, ground) {
    for (i in 1:ncol(deconv)) {
      data = cbind(deconv[, i], ground)
      colnames(data) = c("x", "y")
      mse_val <- mean((data$x - data$y)^2)
      p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_point(color = "blue", size = 0.1, alpha = 0.7) +
        ggplot2::geom_smooth(method = "lm", se = T, color = "red") +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = colnames(ground), y = colnames(deconv)[i]) +
        ggplot2::theme(axis.title.x = ggplot2::element_text(size = 5),
                       axis.text.x = ggplot2::element_text(size = 5),
                       axis.text.y = ggplot2::element_text(size = 5),
                       axis.title.y = ggplot2::element_text(size = 5)) +
        ggplot2::geom_text(ggplot2::aes(x = mean(data$x), y = max(data$y),
                                        label = paste("MSE =", round(mse_val, 4))),
                           size = 2, hjust = 0.5, vjust = -1)
      print(p)
    }
  }
  
  cell_clusters = colnames(groundtruth)
  mse_matrix = data.frame(matrix(ncol = length(deconvolution_combinations),
                                 nrow = length(cell_clusters)))
  rownames(mse_matrix) = cell_clusters
  colnames(mse_matrix) = deconvolution_combinations
  
  cells_discard = c()
  
  for (i in 1:length(cell_clusters)) {
    idx = grep(paste0("_", cell_clusters[i], "$"), colnames(deconvolution))
    if (length(idx) == 0) {
      cells_discard = c(cells_discard, cell_clusters[i])
      next
    }
    deconv = deconvolution[, idx, drop = FALSE]
    ground = groundtruth[, cell_clusters[i], drop = FALSE]
    
    if (scatter && ncol(deconv) != 0) {
      pdf(paste0("Scatter_MSE_", colnames(ground), "_", file_name))
      scatter_plots(deconv, ground)
      dev.off()
    }
    
    for (j in 1:ncol(mse_matrix)) {
      sub_idx = grep(colnames(mse_matrix)[j], colnames(deconv))
      if (length(sub_idx) == 0) {
        mse_matrix[i, j] = NA
      } else {
        mse_matrix[i, j] = (deconv[, sub_idx] - ground)^2
      }
    }
  }
  
  if (length(cells_discard) > 0) {
    mse_matrix = mse_matrix[-which(rownames(mse_matrix) %in% cells_discard), ]
  }
  
  mse_matrix[nrow(mse_matrix) + 1, ] = colMeans(mse_matrix, na.rm = TRUE)
  rownames(mse_matrix)[nrow(mse_matrix)] = "average"
  
  mse_matrix = t(mse_matrix) %>% data.frame() %>%
    dplyr::arrange(dplyr::desc(average)) %>% t() %>% data.frame()
  
  mse_df <- reshape2::melt(mse_matrix)
  mse_df = mse_df %>%
    dplyr::mutate(Cells = rep(rownames(mse_matrix), ncol(mse_matrix)))
  
  g <- mse_df %>% ggplot2::ggplot(ggplot2::aes(Cells, variable,
                                               fill = value, label = round(value, 4))) +
    ggplot2::geom_tile() +
    ggplot2::labs(x = NULL, y = NULL, fill = "Mean Squared Error",
                  title = file_name) +
    ggplot2::scale_fill_gradient(low = "#FCAAAA", high = "white") +
    ggplot2::geom_text(size = 4) +
    ggplot2::theme_classic() +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggpubr::rotate_x_text(angle = 45) +
    ggplot2::theme(axis.text.x = ggtext::element_markdown(),
                   axis.text.y = ggtext::element_markdown())
  
  if (plot) {
    ggplot2::ggsave(
      filename = paste0("Results/Benchmark_MSE_plot_", file_name, ".svg"),
      plot = g,
      width = width,
      height = height,
      device = "svg"
    )
  }
  
  return(mse_matrix)
}
