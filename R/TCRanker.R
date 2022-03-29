# Function TCRanker
#' Rank TCRs by clonotypes for their level in a certain state
#'
#' Calculates clonal scores of all different clonotypes for their level in a
#'  certain state, exhaustion level by default, and provides their final
#'  ranking. The ranking state could be customized by providing corresponding
#'  gene signature. For more details please see the help vignette:
#' \code{vignette("help", package = "TCRanker")}.
#'
#' @param query \code{SingleCellExperiment} object, \code{Seurat} object or
#'  Matrix/DF with query data.
#' @param tcr Vector containing TCRs clonotype info. or a simple entry name
#'  (only when query is of \code{SingleCellExperiment / Saurat} object).
#' @param signature The gene signature that represents a cell status.
#'  By default is the signature for exhaustion level of CD8+ T-cells.
#' @param assay Name of expression data assay, only valid for
#'  \code{SingleCellExperiment / Saurat} object. By default
#'  "counts" for \code{SingleCellExperiment} and
#'  "RNA" for \code{Seurat}.
#' @param plot Logical, to plot the result or not, TRUE by default.
#' @param group Vector or a simple entry name (only when query is of
#'  \code{SingleCellExperiment / Saurat} object). Extra grouping information
#'  to be aggregated and included in the output. Same clonotype in different
#'  groups would be aggregated separately. Even thought optional, it's still
#'  recommended to be included.
#' @param FUN Function used to aggregate scores of the same clonotype. It could
#'  be mean, median or customized functions that take a numeric vector or list
#'  as input. By default, it uses "mean".
#' @param minClonSize Threshold of clonal size that would be taken into account,
#'  5 by default.
#' @param filterCell Name of the sub cell type to filter using \code{scGate},
#'  "CD8" by default. Set to "None" to disable (unfinished function)
#' @param species Charactor, "mouse" or "human", optional.
#' @return a data frame including clonal scores and ranking
#'  (and annotations, if \code{annotation} was included)
#' @examples
#' TCRanker(query, tcr, plot=TRUE)
#' @export TCRanker

TCRanker <- function(query, tcr, signature="default", assay=NULL, plot=TRUE,
                     group="none", FUN="mean", minClonSize=5,
                     filterCell="CD8T", species="auto", ...) {
}

setGeneric("TCRanker")