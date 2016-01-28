#' Plots boxplots of core score distributions by condition and by stability for red and/or blue nodes
#' 
#' This function tests whether core scores of stable nodes are higher than those of variable nodes
#' @param x first output of \code{\link{condor.qscore}}
#' @param y second output of \code{\link{condor.qscore}}
#' @param by character string indicating whether to determine stability by 'row', 'column', or 'both'
#' @param label.x Label of first condition
#' @param label.y Label of second condition
#' @param red.name Name for red nodes
#' @param blue.name Name for blue nodes
#' @param scale.log TRUE/FALSE - If TRUE, plot log core scores
#' @param nsamp Number of permutation tests to run, passed to \code{\link{condor.core.enrich}}
#' @param plot.type character string indicating whether to plot violin plots or points
#' @param node.type character string indicating whether to plot core score comparisons for red nodes, blue nodes, or both
#' @param notch TRUE/FALSE - Add notches to boxplots
#' @return ggplot object
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.compare.qscores <- function(x, y, by=c("column","row","both"),
                                   label.x="Condition x", label.y="Condition y",
                                   red.name="red", blue.name="blue", scale.log=FALSE,
                                   nsamp=1e4, plot.type=c("violin", "points"),
                                   node.type=c("both","red", "blue"), notch=TRUE) {
  by <- match.arg(by)
  plot.type <- match.arg(plot.type)
  node.type <- match.arg(node.type)
  
  ylabel <- ifelse(scale.log, "log(core score)", "Core score")
  qscore.red <- calc.qscore.stability(x, y, type="red", by=by)
  qscore.blue <- calc.qscore.stability(x, y, type="blue", by=by)
  
  perm.test <- function(qscore, type) {
    d <- melt(qscore, measure.vars = c("Q.x", "Q.y"))
    stable.genes <- qscore$names[which(qscore$is.stable)]
    x <- qscore[,c("names","com.x","Q.x")]
    y <- qscore[,c("names","com.y","Q.y")]
    
    x.p <- condor.core.enrich(stable.genes, x, perm=TRUE, nsamp=nsamp)$perm.pvals
    y.p <- condor.core.enrich(stable.genes, y, perm=TRUE, nsamp=nsamp)$perm.pvals
    
    df <- data.frame(variable=c("Q.x", "Q.y"),
                     p=c(x.p$wilcox.perm, y.p$wilcox.perm),
                     ks=c(x.p$ks.perm, y.p$ks.perm))
    d <- merge(d, df, by="variable")
    
    d$p <- paste("Wilcoxon p-value = ", format(d$p, scientific=T, digits=2))
    d$ks <- paste("KS p-value = ", format(d$ks, scientific=T, digits=2))
    d$p <- paste(d$p, d$ks, sep='\n')
    d$type <- ifelse(type=="red", red.name, blue.name)
    return(d)
  }
  
  d.red <- d.blue <- data.frame()
  if (sum(qscore.red$is.stable) < nrow(qscore.red)) {
    if (node.type != "blue") {
      d.red <- perm.test(qscore.red, "red")  
    }
  } else {
    warning("No reds are assigned to different communities between conditions")
  }
  if (sum(qscore.blue$is.stable) < nrow(qscore.blue)) {
    if (node.type != "red") {
      d.blue <- perm.test(qscore.blue, "blue")      
    }
  } else {
    warning("No blues are assigned to different communities between conditions")
  }
  d <- rbind(d.red, d.blue)
  if (nrow(d)==0) {
    stop("No reds nor blues are assigned to different communities between conditions. Exiting.")
  }
  
  if (scale.log) {
    d$value <- log(d$value)
  }
  if (plot.type=="violin") {
    ggplot(d, aes(x=variable, y=value, color=is.stable, fill=is.stable, label=p), environment=environment()) +
      geom_boxplot(width=0.1, notch=notch, outlier.colour=NA, position=position_dodge(width=0.8), fill="white") +
      geom_violin(position=position_dodge(width=0.8), adjust=0.5, alpha=0.3) +
      scale_color_manual(name="", labels=c("Variable", "Stable"), values=c("seagreen4", "gray20"), guide=FALSE) +
      scale_x_discrete(labels=c(label.x, label.y)) +
      scale_fill_manual(name="", labels=c("Variable", "Stable"), values=c("seagreen4", "gray20")) +
      xlab("Condition") + ylab(ylabel) +
      geom_text(data=d[!duplicated(paste0(d$variable, d$type)),], aes(x=variable, y=max(d$value) + 1), vjust=1, position=position_dodge(width=0.9), size=10*5/14) +
      facet_grid(type~., scales="free")
  } else {
    ggplot(d, aes(x=variable, y=value, fill=is.stable, color=is.stable, label=p), environment=environment()) +
      geom_boxplot(fill="white", outlier.colour=NA, notch=notch,
                   position=position_dodge(width=0.9)) +
      geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.3) +
      scale_color_manual(name="", labels=c("Variable", "Stable"), values=c("seagreen4", "gray20")) +
      scale_fill_discrete(guide=FALSE) +
      scale_x_discrete(labels=c(label.x, label.y)) +
      xlab("Condition") + ylab(ylabel) +
      ggtitle("Core score distributions by condition and stability") +
      geom_text(data=d[!duplicated(paste0(d$variable, d$type)),], aes(x=variable, y=max(d$value)), vjust=1, position=position_dodge(width=0.9), size=10*5/14) +
      facet_grid(type~., scales="free")
  }
}