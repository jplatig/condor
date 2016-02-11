## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  message=FALSE,
  warning=FALSE
)

## ------------------------------------------------------------------------
library(condor)
library(igraph)
library(ggplot2)

## ------------------------------------------------------------------------
r = c(1,1,1,2,2,2,3,3,3,4,4);
b = c(1,2,3,1,2,4,2,3,4,3,4);
reds <- c("Alice","Sue","Janine","Mary")
blues <- c("Bob","John","Ed","Hank")
elist <- data.frame(red=reds[r], blue=blues[b])

## ------------------------------------------------------------------------
condor.object <- create.condor.object(elist)

## ------------------------------------------------------------------------
names(condor.object)

## ---- results="hide"-----------------------------------------------------
condor.object <- condor.cluster(condor.object)

## ------------------------------------------------------------------------
print(condor.object$red.memb)
print(condor.object$blue.memb)

## ------------------------------------------------------------------------
gtoy = graph.edgelist(as.matrix(elist),directed=FALSE)
set.graph.attribute(gtoy, "layout", layout.kamada.kawai(gtoy))
V(gtoy)[c(reds,blues)]$color <- c(rep("red",4),rep("blue",4))

## ---- fig.height=4, fig.width=4------------------------------------------
plot(gtoy,vertex.label.dist=2)

## ---- results="hide"-----------------------------------------------------
condor.object <- condor.qscore(condor.object)

## ---- fig.height=4-------------------------------------------------------
q_women <- condor.object$qscores$red.qscore
core_stats <- suppressWarnings(condor.core.enrich(test_nodes=c("Alice","Mary"),
                                                  q=q_women,perm=TRUE,plot.hist=TRUE))

## ---- results="hide"-----------------------------------------------------
data(small1976)
condor.object <- create.condor.object(small1976)
condor.object <- condor.cluster(condor.object, project=FALSE)

