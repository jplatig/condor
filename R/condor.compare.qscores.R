#' Plots boxplots of core score distributions by condition and by stability
#' 
#' This function tests whether core scores of stable nodes are higher than those of variable nodes
#' @param qscore Output of \code{\link{calc.qscore.stability}}
#' @param lab.x Label of first condition
#' @param lab.y Label of second condition
#' @param scale.log TRUE/FALSE - If TRUE, plot log core scores
#' @param nsamp Number of permutation tests to run, passed to \code{\link{condor.core.enrich}}
#' @return ggplot object
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.compare.qscores <- function(qscore, lab.x="Condition x", lab.y="Condition y", scale.log=FALSE, nsamp=1e4) {
  d <- melt(qscore, measure.vars = c("Q.x", "Q.y"))
  stable.genes <- qscore$names[which(qscore$stable)]
  x <- qscore[,c("names","com.x","Q.x")]
  y <- qscore[,c("names","com.y","Q.y")]
  
  x.p <- condor.core.enrich(stable.genes, x, perm=TRUE, nsamp=nsamp)$perm.pvals
  y.p <- condor.core.enrich(stable.genes, y, perm=TRUE, nsamp=nsamp)$perm.pvals
  
  df <- data.frame(variable=c("Q.x", "Q.y"),
                   p=c(x.p$wilcox.perm, y.p$wilcox.perm),
                   ks=c(x.p$ks.perm, y.p$ks.perm))
  d <- merge(d, df, by="variable")
  
  d$p <- paste("Wilcox p:", format(d$p, scientific=T, digits=2))
  d$ks <- paste("KS p:", format(d$ks, scientific=T, digits=2))
  d$p <- paste(d$p, d$ks, sep='\n')
  
  ylabel <- "core score"
  if (scale.log) {
    d$value <- log(d$value)
    ylabel <- "log(core score)"
  }
  
  ggplot(d, aes(variable, value, fill=stable, label=p)) +
    geom_boxplot(notch=TRUE) +
    scale_fill_discrete(guide=guide_legend(title=NULL),
                        labels=c("variable", "stable")) +
    scale_x_discrete(labels=c(lab.x, lab.y)) +
    xlab("condition") + ylab(ylabel) +
    ggtitle("Core score distributions by condition and stability") +
    geom_text(data=d[c(1, nrow(d)),], y=max(d$value), vjust=1)
}