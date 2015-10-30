#' Wrapper for core condor functions
#' 
#' Given an input edge list, this function performs community structure clustering,
#' calculates core scores, and generates the network link heatmap.
#' @param return.gcc if TRUE, returns the giant connected component
#' @param cs.method is a string to specify which unipartite community 
#' structure algorithm should be used for the seed clustering. 
#' Options are \code{LCS} (\code{\link[igraph]{multilevel.community}}), 
#' \code{LEC} (\code{\link[igraph]{leading.eigenvector.community}}), 
#' \code{FG} (\code{\link[igraph]{fastgreedy.community}}).
#' @param project Provides options for initial seeding of the bipartite 
#' modularity maximization.
#' If TRUE, the nodes in the first column of \code{condor.object$edges} 
#' are projected and clustered using \code{cs.method}. If FALSE, the 
#' complete bipartite network is clustered using the unipartite clustering 
#' methods listed in \code{cs.method}.
#' @param ... arguments passed to heatmap.2 in \code{\link{condor.plot.heatmap}}
#' @return \code{condor.object}
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- condor(elist)
#' @import igraph
#' @import Matrix
#' @export
#' 
condor <- function(edgelist,return.gcc=TRUE,cs.method="LCS",project=TRUE, ...){
  condor.object <- create.condor.object(edgelist, return.gcc)
  condor.object <- condor.cluster(condor.object, cs.method, project)
  condor.object <- condor.qscore(condor.object)
  condor.plot.heatmap(condor.object, ...)
  return(condor.object)
}