#' Plots overlap between communities for two condor objects with identical nodes
#' 
#' @param x First condor.object membership data.frame
#' @param y Second condor.object membership data.frame
#' @param main Title of plot
#' @param reorder If FALSE, do not reorder rows and columns
#' @return matrix containing overlap between each pair of communities
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.plot.community.overlap <- function(x, y, main="", xlab="Condition 1", ylab="Condition 2", reorder=TRUE) {
  xs <- split(x, x$com)
  ys <- split(y, y$com)
  z <- matrix(nrow=max(x$com), ncol=max(y$com))
  for (i in 1:nrow(z)) {
    for (j in 1:ncol(z))
      z[i, j] <- length(intersect(xs[[i]]$names, ys[[j]]$names))/length(unique(x$names))
  }
  heatmap.2(z, Rowv=reorder, trace="none", main=main,
            xlab=xlab, ylab=ylab,
            col=colorpanel(10, "white", "black"),
            breaks=sort(c(0.1,seq(0, max(z),length.out=10))))
  return(z)
}
