#' Calculate Neighborhood Center Matrix
#'
#' This function calculates a matrix where each row represents a neighborhood. The columns define
#' principal component coordinates for the neighborhood center. 
#'
#' @param milo_obj A `Milo` class object, containing the neighborhood information, including the neighborhood indices `@nhoodIndex` and the neighborhood connectivity matrix `@nhoods`
#' @param cellPCs A matrix where each row corresponds to a cell and each column corresponds to a principal component. The first 12 PCs are used to evaluate neighborhood centers.
#'
#' @return A matrix where each row represents the mean values of the first 12 PCs for all cells in a specific neighborhood. Row names correspond to neighborhood identifiers, and columsn correspond to  mean PC values. 
#' @export
calculateNeighborhoodCenterMatrix <- function(milo_obj, cellPCs) {
  indices <- 1:length(milo_obj@nhoodIndex)
  names(indices) <- as.character(unlist(milo_obj@nhoodIndex))
  neighborhoods <- milo_obj@nhoods
  NeighborhoodCenters <- lapply(1:ncol(neighborhoods),function(i){
    NN <- neighborhoods[,i]
    NN <- names(NN[NN!=0])
    return(colMeans(cellPCs[rownames(cellPCs) %in% NN,1:12]))
  }) |> do.call(what = "rbind")
  
  rownames(NeighborhoodCenters) <- colnames(neighborhoods)
  return(NeighborhoodCenters)
}