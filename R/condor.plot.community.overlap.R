#' Plots overlap between communities for two condor objects with identical nodes
#' 
#' @param x first output of \code{\link{condor.cluster}},  
#' \code{\link{condor.modularity.max}}, or \code{\link{condor.qscore}}
#' @param y second output of \code{\link{condor.qscore}}, 
#' \code{\link{condor.modularity.max}}, or \code{\link{condor.qscore}}
#' @param type character indicating whether overlaps should be determined for blue or red nodes
#' @param scale character indicating if the values should be scaled proportionally in either the row direction or the column direction, or none. The default is \code{"none"}.
#' @return produces a ggplot object
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.plot.community.overlap <- function(x, y, type=c("red", "blue"), scale=c("none","row","column")) {
  if(is.null(x$red.memb) | is.null(x$blue.memb) | is.null(y$red.memb) | is.null(y$blue.memb)){
    stop("Community Memberships missing. Run condor.cluster or condor.modularity.max first!")
  }
  type <- match.arg(type)
  scale <- match.arg(scale)
  x <- x[[sprintf("%s.memb", type)]]
  y <- y[[sprintf("%s.memb", type)]]
  
  # calculate overlap
  xs <- split(x, x$com)
  ys <- split(y, y$com)
  overlap <- matrix(nrow=max(x$com), ncol=max(y$com))
  overlap.scaled <- overlap
  for (i in 1:nrow(overlap)) {
    for (j in 1:ncol(overlap)) {
      overlap[i, j] <- length(intersect(xs[[i]]$names, ys[[j]]$names))
      overlap.scaled[i, j] <- overlap[i, j] / length(unique(x$names))
    }
  }
  if (scale=="row") {
    overlap.scaled <- overlap.scaled/rowSums(overlap.scaled)
  }
  if (scale=="column") {
    overlap.scaled <- t(t(overlap.scaled)/colSums(overlap.scaled))
  }
  d <- melt(overlap.scaled)
  colnames(d)[1:2] <- c("Var1", "Var2")
  d1 <- melt(overlap)
  d$n <- d1$value
  d$n[d$n==0] <- ""
  p <- ggplot(d, aes(factor(Var2), factor(Var1), fill=value)) +
    geom_tile() +
    # geom_text is in mm = pt * 5/14
    geom_text(aes(label=n), size=10*5/14) +
    xlab("Condition y community") +
    ylab("Condition x community")
  if (scale!="none") {
    p <- p + scale_fill_gradient(name=sprintf("Proportion of\n%s in cell", scale), low="white", high="red")
  } else {
    p <- p + scale_fill_gradient(name="Overlap", low="white", high="red")
  }
  p
}
