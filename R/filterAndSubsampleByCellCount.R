#' Removes individuals from the Seurat object with too few cells, and subsamples individuals with too many cells.
#'
#' Individuals are excluded if their number of cells (cellcount) is in the bottom quartile. Conversely, individuals with
#' more than median + (2*IQR) cells are subsampled. For subsampling, the number of cells to retain is drawn
#' from the distribution of cell counts from the remaining individuals.
#' 
#' @param seurat_obj A Seurat object with a meta.data field called 'projid', corresponding to the individual identifier
#' @return a filtered Seurat object.
#' @export
filterAndSubsampleByCellCount <- function(seurat_obj) {
  cellcounts <- table(seurat_obj@meta.data$projid) |> as.data.frame()
  
  # Plot a histogram of the cell counts (DEBUGGING)
  # hist(cellcounts$Freq, 
  #      breaks = 30,   # You can adjust the number of breaks (bins) as needed
  #      main = "Histogram of Cell Counts per Patient",
  #      xlab = "Number of Cells",
  #      ylab = "Frequency",
  #      col = "lightblue",  # Optional: set the color of the bars
  #      border = "black")   # Optional: set the color of the bar borders
  
  min <- round(summary(cellcounts$Freq)[2]) |> unname() # set min threshold as the first quartile
  max <- round(median(cellcounts$Freq) + IQR(cellcounts$Freq) * 2) # set max threshold as median + 2*IQR
  IDs_too_high <- cellcounts[cellcounts$Freq > max,"Var1"] #IDs of individuals with too many cells
  meta_too_High <- raw_cells@meta.data[raw_cells@meta.data$projid %in% IDs_too_high,] # get the corresponding metadata   
  cellcounts <- cellcounts[cellcounts$Freq >= min & cellcounts$Freq <= max,] # remove unwanted individuals from cellcounts
  metaKeep <- raw_cells@meta.data[raw_cells@meta.data$projid %in% cellcounts$Var1,] # keep metadata from individuals in cellcounts (which is now filtered to the individuals we are keeping)
  
  ##Down sample IDs with too many cells
  meta_too_High <- lapply(unique(meta_too_High$projid),function(projid){ # loop through each unique donor in meta_too_high
    sub_too_high <- meta_too_High[meta_too_High$projid == projid,]
    sub_too_high <- sub_too_high[sample(rownames(sub_too_high),size = sample(cellcounts$Freq,1)),]
    return(sub_too_high)
  }) |> do.call(what = "rbind")
  metaKeep <- rbind(metaKeep,meta_too_High)
  countsKeep <- raw_cells@assays$RNA@counts[,rownames(metaKeep)]
  countsKeep <- countsKeep[rowSums(countsKeep)>0,]
  seurat_obj <- CreateSeuratObject(counts = countsKeep)
  seurat_obj$projid <- metaKeep$projid
  seurat_obj$cell_type <- metaKeep$cell_type_high_resolution
  return(seurat_obj)
}
