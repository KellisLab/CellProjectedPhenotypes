% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateCellScores.R
\name{calculateCellScores}
\alias{calculateCellScores}
\title{Calculate Cell Phenotype Scores}
\usage{
calculateCellScores(
  milo_obj,
  neighborhoodPhenotypeAssociations,
  cellNeighborhoodDistanceMatrix,
  col_name,
  varsToTest
)
}
\arguments{
\item{milo_obj}{A \code{Milo} class object, containing the neighborhood connectivity matrix (\verb{@nhoods}) and cell metadata (\verb{@colData}).}

\item{neighborhoodPhenotypeAssociations}{A list where each element corresponds to a phenotype and contains a data frame with neighborhood-level associations. Produced by \code{calculateNeighborhoodCenterMatrix}.}

\item{cellNeighborhoodDistanceMatrix}{A matrix representing the Euclidean distances between cells and the neighborhoods in which they are included, output of \code{cellNeighborhoodDistanceMatrix}.}

\item{col_name}{The name of the column in \code{milo_obj@colData} that contains the cell identifiers.}

\item{varsToTest}{A character vector of phenotypes for which to calculate cell scores.}
}
\value{
A data frame where each column corresponds to a phenotype, and each row corresponds to a cell, containing the calculated phenotype scores for each cell.
}
\description{
This function calculates phenotype scores for each cell based on the association of its neighborhoods with particular phenotypes.
The scores are derived from the log fold change (logFC) values in the \code{neighborhoodPhenotypeAssociations} data structure and
are adjusted by the distances between cells and neighborhood centers.
}
