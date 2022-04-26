## Function rankClonalScores
#' @importFrom UCell ScoreSignatures_UCell
rankClonalScores <- function(expMat, tcrVec, signature="default", groVec="none",
                             exhaustion=TRUE, proliferation=TRUE,
                             species="auto", FUN="mean", minClonSize=5){

    ## Preparing gene signature, Exhaustion signature as default
    if (all(signature=="default")) {
        if (!(exhaustion|proliferation)) {
            stop("Please provide a customized gene signature, ",
                 "or enable at least one default gene signature")
        }
        features <- list()
    } else {
        if (class(signature)=='character') {
            features <- list()
            features$user <- signature
        } else if (class(signature)=='list') {
            if (is.null(names(signature))) {
                stop("signature list should contain names")
            } else features <- signature
        } else stop("Unrecognisable gene signature format")
    }
    
    if (species=="auto") species <- getSpecies(genes=rownames(expMat))
    
    if (exhaustion) {
        if (species == "mouse") {
            message("Using default mouse gene signature.")
            features$exhaustion <- c('Lag3', 'Prf1', 'Havcr2', 'Gzmb', 'Nkg7',
                                     'Ctsd', 'Id2', 'Klrd1', 'Cst7', 'Pdcd1',
                                     'Tnfrsf9', 'Tigit', 'Ctsw', 'Ccl4', 'Cd63',
                                     'Ccl3', 'Ifng', 'Cxcr6', 'Fasl', 'Rbpj',
                                     'Chst12', 'Fam3c', 'Csf1')
        } else if (species == "human") {
            message("Using default human gene signature.")
            features$exhaustion <- c('LAG3', 'PRF1', 'HAVCR2', 'GZMB', 'NKG7', 
                                     'CTSD', 'ID2', 'KLRD1', 'CST7', 'PDCD1',
                                     'TNFRSF9', 'TIGIT', 'CTSW', 'CCL4', 'CD63',
                                     'CCL3', 'IFNG', 'CXCR6', 'FASLG', 'RBPJ',
                                     'CHST12', 'FAM3C', 'CSF1')
        } else stop("Only support human or mouse default gene signature.")
    }
    if (proliferation) {
        if (species == "mouse") {
            features$proliferation <- 
                c(scGate::genes.blacklist.default$Mm$cellCycle.G1S,
                  scGate::genes.blacklist.default$Mm$cellCycle.G2M)
        } else if (species == "human") {
            features$proliferation <- 
                c(scGate::genes.blacklist.default$Hs$cellCycle.G1S,
                  scGate::genes.blacklist.default$Hs$cellCycle.G2M)
        } else stop("Only support human or mouse default gene signature.")
    }
    
    ## Calculating individual scores
    message("Calculating individual scores with UCell")
    individualScores <- ScoreSignatures_UCell(expMat, features=features)
    
    ## Creating full data frame
    fullData <- data.frame(tcrVec, groVec, individualScores)
    colnames(fullData) <- c("clonotype", "group", paste0(names(features),'.score'))
    
    ## Gathering necessary info.
    scores <- aggregate(.~clonotype+group, fullData, FUN=match.fun(FUN)) # clonal score
    sizes <- aggregate(.~clonotype+group, fullData, FUN=length) # clonal size
    colnames(sizes)[3] <- "size"
    sizes <- sizes[sizes$size >= minClonSize,1:3] # screening with minClonSize requirement
    groupSize <- aggregate(size~group,sizes,FUN=sum) # group size
    colnames(groupSize)[2] <- "freq" # storing group size in freq column
    sizes <- merge(sizes, groupSize, by="group") # cbind size and freq
    sizes$freq <- sizes$size/sizes$freq # change freq to percentage
    TCRanking <- merge(sizes, scores, by=c("clonotype","group")) # cbind current table with score
    
    ## Rankings according to each signature score
    rankings <- vapply(X=-TCRanking[,5:ncol(TCRanking)],
                       FUN = rank,
                       FUN.VALUE = numeric(nrow(TCRanking)))
    colnames(rankings) <- paste0(names(features),'.ranking')
    
    TCRanking <- cbind(TCRanking, rankings) # combine ranking to TCRanking
    
    ## Reorder columns
    colOrder <- c(5:(4+ncol(individualScores)))
    colOrder <- c(rbind(colOrder,colOrder+ncol(individualScores)))
    TCRanking <- TCRanking[,c(1:4,colOrder)]
    
    ## remove group column if it wasn't assigned at first
    if (all(TCRanking$group=="none")) {
        TCRanking <- TCRanking[,-2]
    }
    return(TCRanking)
}

## Add aggregating function here below in the future.

## Taking only top "percentage" scores into account | at least "min" cells
top <- function(scores, percentage=0.1, min=5){
    numTop <- round(length(scores)*percentage)
    if(numTop<min) numTop <- min
    score <- mean(tail(sort(scores), numTop))
    return(score)
}

## Taking the 90 percentile as the clonal score
ninetyPT <- function(scores, percentile=0.9){
    return(quantile(scores, probs=percentile))
}
