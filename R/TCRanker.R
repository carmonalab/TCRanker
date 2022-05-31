# Function TCRanker
#' Rank TCRs by clonotypes for their level in a certain state
#'
#' Calculates clonal scores of all different clonotypes for their level in a
#'  certain state, exhaustion level by default, and provides their final
#'  ranking. The ranking state could be customized by providing corresponding
#'  gene signature. For more details please see the help vignette:
#' \code{vignette("help", package = "TCRanker")}.
#'
#' @param query \code{SingleCellExperiment} object or \code{Seurat} object.
#' @param tcr A colname of \code{query} of clonotypes or a vector containing
#'  TCRs clonotype info.
#' @param signature The gene signature that represents a cell status.
#'  By default are the signatures of exhaustion and proliferation level of 
#'  CD8+ T cells. For custom gene signatures, please refer
#'   \href{https://carmonalab.github.io/TCRanker.demo/demo.html}{TCRanker.demo}
#' @param assay Name of expression data assay, By default "counts" for
#'  \code{SingleCellExperiment} and "RNA" for \code{Seurat}.
#' @param group A colname of \code{query} or a vector. Extra group information
#'  to be aggregated and included in the output. Same clonotype in different
#'  groups would be aggregated separately. Even thought optional, it's still
#'  recommended to be included.
#' @param exhaustion Logical, to include exhaustion scores and ranking in the
#'  output. \code{TRUE} by default.
#' @param proliferation Logical, to include proliferation scores and ranking in
#'  the output. \code{TRUE} by default.
#' @param species Charactor, "mouse"/"human", will be auto-detected if omitted.
#' @param FUN Function used to aggregate scores of the same clonotype. It could
#'  be mean, median or customized functions that take a numeric vector or list
#'  as input and return a single numeric. By default, it uses "mean".
#' @param minClonSize Threshold of clonal size that would be taken into account,
#'  5 by default.
#' @param filterCell Name of the sub cell type to filter using \code{scGate},
#'  "CD8T" by default. Set to "none" to disable.
#' @param strictFilter Logical, to exclude impure cells in pure clonotype or
#'  not. Only valid when \code{filterCell} is on. \code{TRUE} by default.
#' @param keepObject Logical, to return the \code{SingleCellExperiment / Saurat}
#'  object after cell filtering or not (returned together with the ranking 
#'  data frame in a list). \code{FALSE} by default.
#' @return a data frame including clonal scores and ranking
#'  (and groups, if \code{group} was specified)
#' 
#' @export TCRanker

TCRanker <- function(query, tcr, signature="default", assay=NULL, group="none", 
                     exhaustion=TRUE, proliferation=TRUE, species="auto", 
                     FUN="mean", minClonSize=5, filterCell="CD8T",
                     strictFilter=TRUE, keepObject=FALSE, ...) {
}

setGeneric("TCRanker")
