#' Classifies each node as either stable or variable between two conditions
#' 
#' Classifies each node as either stable or variable based on the output of \code{\link{condor.plot.community.overlap}}
#' 
#' @param res Output of \code{\link{condor.plot.community.overlap}}
#' @param by character string indicating whether to determine stability by 'row', 'column', or 'both'
#' @return A data frame with modularity of both conditions and stability
#' @export
#' 
calc.qscore.stability <- function(res, by=c("both","row","column")) {
  if (is.null(res$x$Q) | is.null(res$y$Q)) {
    stop("condor.plot.community.overlap must be run on qscores to use this function.")
  }
  by <- match.arg(by)
  qscore <- merge(res$x, res$y, by="names")
  com.map.x <- apply(res$z, 1, function(x){which(x==max(x))})
  com.map.y <- apply(res$z, 2, function(x){which(x==max(x))})
  # get all stable genes
  stable.x <- stable.y <- stable.genes <- c()
  for (i in 1:nrow(res$z)) {
    stable.x1 <- subset(res$x, com==i)$names
    stable.x2 <- subset(res$y, com%in%com.map.x[[i]])$names
    stable.x <- c(stable.x, intersect(stable.x1, stable.x2))
  }
  for (i in 1:ncol(res$z)) {
    stable.y1 <- subset(res$y, com==i)$names
    stable.y2 <- subset(res$x, com%in%com.map.y[[i]])$names
    stable.y <- c(stable.y, intersect(stable.y1, stable.y2))
  }
  if (by=="row") {
    stable.genes <- stable.x
  }
  if (by=="column") {
    stable.genes <- stable.y
  }
  if (by=="both") {
    stable.genes <- intersect(stable.x, stable.y)
  }
  qscore$stable <- with(qscore, names %in% stable.genes)
  return(qscore)
}