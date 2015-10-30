#' Plots overlap between communities for two condor objects with identical nodes
#' 
#' @param x First condor.object membership or qscores data.frame
#' @param y Second condor.object membership or qscores data.frame
#' @param scale character indicating if the values should be scaled in either the row direction or the column direction, or none. The default is \code{"none"}.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... other arguments passed to heatmap.2
#' @return matrix containing overlap between each pair of communities
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.plot.community.overlap <- function(x, y, xlab="Condition y", ylab="Condition x", scale.condition=c("none","row","column"), ...) {
  scale.condition <- match.arg(scale.condition)
  xs <- split(x, x$com)
  ys <- split(y, y$com)
  z <- matrix(nrow=max(x$com), ncol=max(y$com))
  zplot <- z
  for (i in 1:nrow(zplot)) {
    for (j in 1:ncol(zplot)) {
      z[i, j] <- length(intersect(xs[[i]]$names, ys[[j]]$names))
      zplot[i, j] <- z[i, j] / length(unique(x$names))
    }
  }
  if (scale.condition=="row") {
    zplot <- zplot/rowSums(zplot)
  }
  if (scale.condition=="column") {
    zplot <- t(t(zplot)/colSums(zplot))
  }
  heatmap.2(zplot, trace="none", xlab=xlab, ylab=ylab, scale="none",
            col=colorpanel(10, "white", "black"),
            breaks=sort(c(0.1,seq(0, max(zplot),length.out=10))),
            ...)
  res <- list(x=x, y=y, z=z)
  return(res)
}
