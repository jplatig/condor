#' Plots boxplots of core score distributions by condition and by stability
#' 
#' This function tests whether core scores of stable nodes are higher than those of variable nodes
#' @param x First condor.object qscores data.frame
#' @param y Second condor.object qscores data.frame
#' @param z Output of \code{\link{condor.plot.community.overlap}}
#' @param lab.x Label of first condition
#' @param lab.y Label of second condition
#' @return ggplot object
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.compare.qscores <- function(x, y, z, lab.x="Condition 1", lab.y="Condition 2") {
  qscore <- merge(x, y, by="names")
  com.map <- apply(r, 1, function(x){which(x==max(x))})
  # get all stable genes
  stable.genes <- c()
  for (i in 1:nrow(r)) {
    stable.x <- subset(x, com==i)$names
    stable.y <- subset(y, com%in%com.map[[i]])$names
    stable.genes <- c(stable.genes, intersect(stable.x, stable.y))
  }
  qscore$stable <- with(qscore, names %in% stable.genes)
  d <- melt(qscore, measure.vars = c("Q.x", "Q.y"))
  
  x.p <- wilcox.test(value ~ stable, data=subset(d, variable=="Q.x"))$p.value
  y.p <- wilcox.test(value ~ stable, data=subset(d, variable=="Q.y"))$p.value
  
  x.ks <- ks.test(subset(qscore, stable==T)$Q.x, subset(qscore, stable==F)$Q.x, alternative="less")$p.value
  y.ks <- ks.test(subset(qscore, stable==T)$Q.y, subset(qscore, stable==F)$Q.y, alternative="less")$p.value
  
  df <- data.frame(variable=c("Q.x", "Q.y"), p=c(x.p, y.p), ks=c(x.ks, y.ks))
  d <- merge(d, df, by="variable")
  
  d$p <- paste("Wilcox p:", format(d$p, scientific=T, digits=2))
  d$ks <- paste("KS p:", format(d$ks, scientific=T, digits=2))
  d$p <- paste(d$p, d$ks, sep='\n')
  
  ggplot(d, aes(variable, value, fill=stable, label=p)) + geom_boxplot() +
    scale_fill_discrete(guide=guide_legend(title=NULL),
                        labels=c("variable", "stable")) +
    scale_x_discrete(labels=c(lab.x, lab.y)) +
    xlab("condition") + ylab("core score") +
    ggtitle("Core score distributions by condition and stability") +
    geom_text(data=d[c(1, nrow(d)),], y=max(d$value), vjust=1)
}