## Function rankClonalScores
#' @importFrom UCell ScoreSignatures_UCell
rankClonalScores <- function(expMat, tcrVec, signature="default", groVec="none",
                             species="auto", FUN="mean", minClonSize=5){

    ## Preparing gene signature, Exhaustion signature as default
    if (all(signature=="default")) {
        if (species=="auto") {
            species <- getSpecies(genes=rownames(expMat))
        }
        if (species == "mouse"){
            message("Using default mouse gene signature.")
            signature <- c('Lag3', 'Prf1', 'Havcr2', 'Gzmb', 'Nkg7','Ctsd',
                           'Klrd1', 'Id2', 'Cst7', 'Pdcd1', 'Tnfrsf9', 'Tigit',
                           'Ctsw', 'Ccl4', 'Cd63', 'Ccl3', 'Ifng', 'Cxcr6',
                           'Fasl', 'Rbpj', 'Chst12', 'Fam3c', 'Csf1')
        } else if (species == "human") {
            message("Using default human gene signature.")
            signature <- c('LAG3', 'PRF1', 'HAVCR2', 'GZMB', 'NKG7', 'CTSD',
                           'KLRD1', 'ID2', 'CST7', 'PDCD1', 'TNFRSF9', 'TIGIT',
                           'CTSW', 'CCL4', 'CD63', 'CCL3', 'IFNG', 'CXCR6',
                           'FASLG', 'RBPJ', 'CHST12', 'FAM3C', 'CSF1')
        } else message("Only support human or mouse gene signature.")
    }

    ## Calculating scores
    message("Calculating individual scores with UCell")
    individualScores <- ScoreSignatures_UCell(expMat,
                                              features=list(Tex=signature))
    ## Creating full data frame
    fullData <- data.frame(tcrVec, individualScores, groVec)
    colnames(fullData) <- c("clonotype", "score", "group")
    ## Gathering necessary info.
    scores <- aggregate(score~., fullData, FUN=match.fun(FUN))
    sizes <- aggregate(score~., fullData, FUN=length)
    colnames(sizes)[3] <- "size"
    sizes <- sizes[sizes$size >= minClonSize,]
    groupSize <- aggregate(size~group,sizes,FUN=sum)
    colnames(groupSize)[2] <- "freq"
    sizes <- merge(sizes, groupSize, by="group")
    sizes$freq <- sizes$size/sizes$freq
    TCRanking <- merge(scores, sizes, by=c("clonotype","group"))
    TCRanking$ranking <- rank(-TCRanking$score)
    TCRanking <- TCRanking[order(TCRanking$ranking),] # re-order
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
