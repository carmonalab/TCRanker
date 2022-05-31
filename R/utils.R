# Internal function to get the TCR clonotype vector
getTCR <- function(tcr, expMat) {
    ids <- names(tcr)
    if (is.null(ids)) {
        if (length(tcr) == dim(expMat)[2]) {
            names(tcr) <- colnames(expMat)
            tcrVec <- tcr
        }
        else stop("tcr should either has the same length as the cell number ",
                  "or contain cell IDs as names")
    } else {
        inter <- intersect(ids, colnames(expMat))
        if (length(inter) == 0) {
            stop("None of the tcr names were found in the query.")
        }
        else {
            tcrVec <- tcr[inter]
            if (length(inter) != length(ids)) {
                warning("Some tcr names were not found in the query.")
            }
        }
    }
    return(tcrVec)
}

# Internal function to get group vector
getGroup <- function(group, expMat) {
    if (any(group!="none")) {
        ids <- names(group)
        if (is.null(ids)) {
            if (length(group) == dim(expMat)[2]) {
                names(group) <- colnames(expMat)
                groVec <- group
            }
            else stop("group should either has the same length as the cell ",
                      "number or contain cell IDs as names")
        } else {
            inter <- intersect(ids, colnames(expMat))
            if (length(inter) == 0) {
                stop("None of the group names were found in the query.")
            }
            else {
                groVec <- group[inter]
                if (length(inter) != length(ids)) {
                    warning("Some group names were not found in the query.")
                }
            }
        }
    } else groVec <- rep("none", dim(expMat)[2])

    return(groVec)
}

# Internal function to filter target functional cluster by scGate
#' @import scGate
filterCells <- function(query, tcr, species="mouse", cluster="CD8T", 
                        assay=NULL, strictFilter=TRUE){
    # counts of cells
    ncells <- ncol(query)
    
    if(length(tcr) == 1) {
        query$clonotype <- FetchData(query, tcr) # tcr as colname of meta.data
    } else {
        query$clonotype <- tcr # tcr as vector
    }
    
    # remove cells with NA colotype
    query <- subset(query, cells = Cells(query)[!is.na(query$clonotype)])
    nNAcells <- ncells - ncol(query) # count of cells with NA clonotype 
    # fetch gating model
    models<-suppressMessages(get_scGateDB())
    model <- tryCatch(expr=models[[species]]$generic[[cluster]],
                      error=function(cond)stop("Functional cluster ", cluster,
                                               " not found in ", species,
                                               " models in scGate.",
                                               "Please choose from ",
                                               names(models$species$generic)))
    # scGate
    query <- suppressWarnings(scGate(data=query, model = model,
                                            verbose=FALSE, assay=assay))
    # select pure clonotype (>50% pure cells)
    pure.TCRs <- tapply(query$is.pure=="Pure", query$clonotype, mean) > 0.5
    pure.TCRs <- names(subset(pure.TCRs,subset=pure.TCRs))
    
    if(strictFilter) { # exclude impure cells in pure clonotype
        query <- subset(query, 
                        subset=(clonotype %in% pure.TCRs & is.pure=="Pure"))
    } else { # include impure cells in pure clonotype
        query <- subset(query, subset=(clonotype %in% pure.TCRs))
    }
    
    # counts of remaining cell
    ncells.keep <- ncol(query)
    if (ncells.keep > 0) {
        message(ncells - ncells.keep, " (", nNAcells, " with NA clonotype, ",
                ncells - ncells.keep - nNAcells, " impure) out of ", ncells,
                " ( ", round(100*(ncells-ncells.keep)/ncells,digits = 2),
                "% ) cells removed.")
    } else {
        stop("Stopping. All cells were removed by cell filter!")
    }
    # clear up meta.data
    to_remove <- grep("is.pure", colnames(query@meta.data))
    to_remove <- c(to_remove, 
                   grep("_UCell$", colnames(query@meta.data), perl=T))

    query@meta.data <- query@meta.data[,-to_remove]
    return(query)
}

# Internal function to determine the species
getSpecies <- function(genes, table=Hs2Mm.convert.table) {

    g.mm <- length(intersect(genes, table$Gene.MM))
    g.hs1 <- length(intersect(genes, table$Gene.stable.ID.HS))
    g.hs2 <- length(intersect(genes, table$Gene.HS))
    gg <- c(g.mm, g.hs1, g.hs2)

    if (max(gg)==g.mm) {
        species='mouse'
    } else {
        species='human'
    }
    return(species)
}
