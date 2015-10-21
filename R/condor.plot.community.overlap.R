#' Plots overlap between communities for two condor objects with identical nodes
#' 
#' @param x First condor.object membership data.frame
#' @param y Second condor.object membership data.frame
#' @param main Title of plot
#' @param reorder If FALSE, do not reorder rows and columns
#' @param scale character indicating if the values should be scaled in either the row direction or the column direction, or none. The default is \code{"none"}.
#' @return matrix containing overlap between each pair of communities
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.plot.community.overlap <- function(x, y, main="", xlab="Condition y", ylab="Condition x", reorder=TRUE, scale="none") {
  xs <- split(x, x$com)
  ys <- split(y, y$com)
  z <- matrix(nrow=max(x$com), ncol=max(y$com))
  for (i in 1:nrow(z)) {
    for (j in 1:ncol(z)) {
      z[i, j] <- length(intersect(xs[[i]]$names, ys[[j]]$names))/length(unique(x$names))
    }
  }
  if (scale=="row") {
    z <- z/rowSums(z)
  }
  if (scale=="column") {
    z <- t(t(z)/colSums(z))
  }
  heatmap.2(z, Rowv=reorder, trace="none", main=main,
            xlab=xlab, ylab=ylab, scale="none",
            col=colorpanel(10, "white", "black"),
            breaks=sort(c(0.1,seq(0, max(z),length.out=10))))
  return(z)
}
