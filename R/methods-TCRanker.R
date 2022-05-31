## TCRanker method for unsupported query class
TCRanker_ANY <- function(query) {
    stop("Unrecognized input query format.",
         "Please use Seurat or SCE object as the input format of query;")
}

## TCRanker method for SCE object query
#' @importFrom SummarizedExperiment assay colData 
#' @importFrom SingleCellExperiment SingleCellExperiment
TCRanker_SCE <- function(query, tcr, signature="default", assay="counts",
                         group="none", exhaustion=TRUE, proliferation=TRUE,
                         species="auto", FUN="mean", minClonSize=5,
                         filterCell="CD8T", strictFilter=TRUE, 
                         keepObject=FALSE) {
    ## Filter functional cluster
    if(filterCell != "none"){
        message("Filtering ", filterCell, " Cells. ",
                "Set filterCell='none' to disable pre-filtering")
        if(species=="auto") species <- getSpecies(genes=rownames(query))
        query <- as.Seurat(query, counts = assay)
        query <- filterCells(query, tcr, species=species, cluster=filterCell,
                             assay="RNA", strictFilter)
        query <- as.SingleCellExperiment(query)
    }

    ## Get expression matrix
    if(is.null(assay)){ # unspecified assay
        stop("Please specify an assay.")
    } else { # specified assay
        expMat <- tryCatch(expr=assay(query,assay),
                           error=function(cond)stop("Assay ", assay,
                                                    " not found in query"))
        message("Using assay ", assay)
    }

    ## Get TCR vector
    if(length(tcr)==1) {
        tcrVec <- tryCatch(expr=colData(query)[tcr],
                           error=function(cond)stop("TCR colData ", tcr,
                                                    " not found in query"))
        message("Using TCR colData ", tcr)
    } else tcrVec <- getTCR(tcr, expMat)

    ## Get group vector
    if(length(group)==1 & group!="none") {
        groVec <- tryCatch(expr=colData(query)[group],
                            error=function(cond)stop("Meta.data ", group,
                                                     " not found in query"))
        message("Using group meta.data ", group)
    } else groVec <- getGroup(group, expMat)

    TCRanking <- rankClonalScores(expMat, tcrVec, signature, groVec, exhaustion,
                                  proliferation, species, FUN, minClonSize)
    
    if(keepObject){
        result_list <- list("TCRanking" = TCRanking, "filteredQuery" = query)
        return(result_list) 
    }

    message("Done")
    return(TCRanking)
}

## TCRanker method for Seurat object query
#' @importFrom SeuratObject GetAssayData FetchData
TCRanker_Seurat <- function(query, tcr, signature="default", assay="RNA",
                            group="none", exhaustion=TRUE, proliferation=TRUE,
                            species="auto", FUN="mean", minClonSize=5,
                            filterCell="CD8T", strictFilter=TRUE, 
                            keepObject=FALSE) {
    ## Filter functional cluster
    if(filterCell != "none"){
        message("Filtering ", filterCell, " Cells. ",
                "Set filterCell='none' to disable pre-filtering")
        if(species=="auto") species <- getSpecies(genes=rownames(query))
        query <- filterCells(query, tcr, species=species, cluster=filterCell,
                             assay=assay, strictFilter)
    }

    ## Get expression matrix
    if(is.null(assay)){ # unspecified assay
        stop("Please specify an assay")
    } else { # specified assay
        expMat <- tryCatch(expr=GetAssayData(query, assay=assay),
                           error=function(cond)stop("Assay ", assay,
                                                    " not found in query"))
        message("Using assay ", assay)
    }

    ## Get TCR vector
    if(length(tcr) == 1) {
        tcrVec <- tryCatch(expr=FetchData(query, tcr),
                           error=function(cond)stop("Meta.data ", tcr,
                                                    " not found in query"))
        message("Using TCR meta.data ", tcr)
    } else tcrVec <- getTCR(tcr, expMat)


    ## Get group vector
    if(length(group)==1 & group!="none") {
        groVec <- tryCatch(expr=FetchData(query, group),
                            error=function(cond)stop("Meta.data ",
                                                     group,
                                                     " not found in query"))
        message("Using group meta.data ",group)
    } else groVec <- getGroup(group, expMat)

    TCRanking <- rankClonalScores(expMat, tcrVec, signature, groVec, exhaustion,
                                  proliferation, species, FUN, minClonSize)
    
    if(keepObject){
        result_list <- list("TCRanking" = TCRanking, "filteredQuery" = query)
        return(result_list) 
    }
    
    message("Done")
    return(TCRanking)
}

setMethod("TCRanker", signature(query="ANY"), TCRanker_ANY)
setMethod("TCRanker", signature(query="SingleCellExperiment"), TCRanker_SCE)
setMethod("TCRanker", signature(query="Seurat"), TCRanker_Seurat)
