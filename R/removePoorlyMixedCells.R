#' Remove cells with poor mixture score from Seurat object
#' 
#' This function calculates the Local Inverse Simpson's Index (LISI)
#' as a metric of cell mixture. Specifically, LISI measures the effective sample
#' size in the local embedding space. Poorly mixed space is defined by the 
#' threshold parameter (recommended value at 3). 
#' 
#' Note that this function contains a commented debug code block for visualizing 
#' the distribution of LISI scores across cells.
#'
#' @param seurat_obj A Seurat object for which you want poorly mixed cells removed.
#' @param cellPCs A matrix of cells (rows) and PCA embeddings (columns). Can be extracted with seurat_obj@reductions$pca@cell.embeddings[,1:12]
#' @param col_name Define the column name to evaluate mixture. Can be individual, or condition, or a covariate.
#' @param threshold A positive real number, defining the minimum LISI required for cell retention. 
#'
#' @importFrom lisi compute_lisi
#' @return A filtered Seurat object
#' @export
removePoorlyMixedCells <- function(seurat_obj, cellPCs, col_name, threshold) {
  lisi_result <- compute_lisi(cellPCs, seurat_obj@meta.data, col_name)
  lisi_result <- as.data.frame(lisi_result)
  lisi_result$ids <- rownames(lisi_result)
  toRemove <- lisi_result[lisi_result[[col_name]] < threshold,"ids"]
  seurat_obj <- seurat_obj[,!colnames(seurat_obj) %in% toRemove]
  # DEBUG: print distribution of LISI scores across cells
  # hist(lisi_res$projid,
  #            breaks = 30,
  #            main = "Histogram of LISI per cell",
  #            xlab = "LISI",
  #            ylab = "Frequency",
  #            col = "lightblue",  # Optional: set the color of the bars
  #            border = "black")   # Optional: set the color of the bar borders
  return(seurat_obj)
}
