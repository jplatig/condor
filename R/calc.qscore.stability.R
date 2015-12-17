#' Classifies each node as either stable or variable between two conditions
#' 
#' Classifies each node as either stable (remaining in the most shared
#' community between conditions) or variable (switching communities)
#' 
#' @param x first output of \code{\link{condor.qscore}}
#' @param y second output of \code{\link{condor.qscore}}
#' @param type character indicating whether overlaps should be determined for blue or red nodes
#' @param by character indicating whether to determine stability by 'row', 'column', or 'both'.
#' Default is "column"
#' @return data.frame with core scores of both conditions and stability
#' @export
#' 
calc.qscore.stability <- function(x, y, type=c("red", "blue"), by=c("column","row","both")) {
  if (is.null(x$qscores) | is.null(y$qscores)) {
    stop("condor.qscore must be run first to use this function.")
  }
  type <- match.arg(type)
  by <- match.arg(by)
  x.memb <- x[[sprintf("%s.memb", type)]]
  y.memb <- y[[sprintf("%s.memb", type)]]
  # calculate overlap
  xs <- split(x.memb, x.memb$com)
  ys <- split(y.memb, y.memb$com)
  overlap <- matrix(nrow=max(x.memb$com), ncol=max(y.memb$com))
  for (i in 1:nrow(overlap)) {
    for (j in 1:ncol(overlap)) {
      overlap[i, j] <- length(intersect(xs[[i]]$names, ys[[j]]$names))
    }
  }
  x.qscore <- x$qscores[[sprintf("%s.qscore", type)]]
  y.qscore <- y$qscores[[sprintf("%s.qscore", type)]]
  qscore <- merge(x.qscore, y.qscore, by="names")
  com.map.x <- apply(overlap, 1, function(x){which(x==max(x))})
  com.map.y <- apply(overlap, 2, function(x){which(x==max(x))})
  # get all stable genes
  stable.x <- stable.y <- stable.genes <- c()
  for (i in 1:nrow(overlap)) {
    stable.x1 <- subset(x.qscore, com==i)$names
    stable.x2 <- subset(y.qscore, com%in%com.map.x[[i]])$names
    stable.x <- c(stable.x, intersect(stable.x1, stable.x2))
  }
  for (i in 1:ncol(overlap)) {
    stable.y1 <- subset(y.qscore, com==i)$names
    stable.y2 <- subset(x.qscore, com%in%com.map.y[[i]])$names
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
  qscore$is.stable <- with(qscore, names %in% stable.genes)
  return(qscore)
}