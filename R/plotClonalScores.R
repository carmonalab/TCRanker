
# Function plotClonalScores
#' Show the distribution of clonal scores
#'
#' Plot the clonal ranking of all clonotypes against their size,
#' including their clonal scores.
#'
#' @param TCRanking Dataframe including clonal scores and ranking.
#'  Could be obtained from \code{TCRanker()}
#' @param title title of the plot
#' @param palette Vector of hex color code, optional.
#' @param labelSize Numeric, font size of clonal scores
#' @param textSize Numeric, font size of clonotype sequence
#' @return Plot the clonal ranking of all clonotypes against their size
#' @examples
#' plotClonalScores(TCRanking)
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export plotClonalScores

plotClonalScores <- function(TCRanking, title=NULL, subtitle=NULL, palette=NULL,
                             labelSize=2, textSize=7){
    ## Checking input format
    if (ncol(TCRanking)<4 |
        all(colnames(TCRanking)[1:4] !=
            c("clonotype", "score", 'size', "ranking"))) {
        stop("Wrong input column number or names ")
    }

    ## Reshaping to long format if multiple groups in one clonotype
    if (ncol(TCRanking)>5) {
        TCRanking <- melt(TCRanking,
                          id.vars=colnames(TCRanking)[1:4],
                          measure.vars=colnames(TCRanking)[-(1:4)],
                          variable.name="group", value.name="subsize")
    }

    ## Setting up palette
    if(is.null(palette)){
        palette <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB",
                     "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
    }
    numGroup <- length(unique(TCRanking$group))
    if (numGroup < length(palette)) {
        palette <- palette[1:numGroup]
    } else palette <- rainbow(n=numGroup)

    ## Initiating plot
    g <- ggplot(TCRanking, aes(x=ranking))

    ## Filling the core mapping
    g <- g + aes(y=size)+
        geom_bar(stat="identity")  +
        geom_text(aes(label=round(score, 4)), size=labelSize,
                  hjust='left', nudge_y=1)
    if (ncol(TCRanking)==5) { # if groups exist
        g <- g +
            aes(fill=group) +
            scale_fill_manual(values=palette) +
            theme(legend.position="bottom")
    }

    ## Filling the rest
    if (is.null(title)) title <- "TCR ranking"
    g <- g +
        coord_flip()+theme_light() +
        scale_x_continuous(trans="reverse", position="top",
                           breaks=seq(1, nrow(TCRanking), 5),
                           expand=c(0, 1),
                           sec.axis=dup_axis(
                               breaks=seq(1, nrow(TCRanking), 1),
                               labels=TCRanking[order(TCRanking$ranking),1],
                               name="Clonotype")) +
        scale_y_continuous(position="right", limits=c(0, NA),
                           breaks=seq(0, max(TCRanking$size), 5),
                           expand=expansion(mult=c(0, 0.05)))+
        labs(title=title,
             subtitle=subtitle,
             x="Clonal Ranking(Scores)",
             y="Clonal Size") +
        theme(plot.title=element_text(hjust=0.5),
              plot.subtitle=element_text(hjust=0.5),
              axis.text.y=element_text(size=textSize))
    return(g)
}


