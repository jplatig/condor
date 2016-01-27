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

## ---- include=FALSE------------------------------------------------------
condor.object <- condor.cluster(condor.object)

## ------------------------------------------------------------------------
print(condor.object$red.memb)
print(condor.object$blue.memb)

## ------------------------------------------------------------------------
gtoy = graph.edgelist(as.matrix(elist),directed=FALSE)
set.graph.attribute(gtoy, "layout", layout.kamada.kawai(gtoy))
V(gtoy)[c(reds,blues)]$color <- c(rep("red",4),rep("blue",4))

## ------------------------------------------------------------------------
plot(gtoy,vertex.label.dist=2)

## ----include=FALSE-------------------------------------------------------
condor.object <- condor.qscore(condor.object)

## ------------------------------------------------------------------------
q_women <- condor.object$qscores$red.qscore
core_stats <- suppressWarnings(condor.core.enrich(test_nodes=c("Alice","Mary"),
                                                  q=q_women,perm=TRUE,plot.hist=TRUE))

## ---- include=FALSE------------------------------------------------------
data(small1976)
condor.object <- create.condor.object(small1976)
condor.object <- condor.cluster(condor.object, project=FALSE)

## ------------------------------------------------------------------------
condor.plot.heatmap(condor.object, xlab="Plants", ylab="Pollinators", add.color=TRUE)

## ---- include=FALSE------------------------------------------------------
set.seed(1)
small1976.noisy <- small1976
small1976.noisy[, 3] <- small1976[, 3] + floor(runif(nrow(small1976), -5, 5))
small1976.noisy[which(small1976.noisy[, 3] < 0), 3] <- 0
condor.object.noisy <- create.condor.object(small1976.noisy)
condor.object.noisy <- condor.cluster(condor.object.noisy, project=FALSE)

## ------------------------------------------------------------------------
condor.plot.heatmap(condor.object.noisy, xlab="Plants", ylab="Pollinators", add.color=TRUE)

## ------------------------------------------------------------------------
condor.plot.community.overlap(condor.object, condor.object.noisy) +
  ylab("Original data-set community") + xlab("Noise-added data-set community")

## ----include=FALSE-------------------------------------------------------
condor.object <- condor.qscore(condor.object)
condor.object.noisy <- condor.qscore(condor.object.noisy)

## ------------------------------------------------------------------------
condor.compare.qscores(condor.object, condor.object.noisy, nsamp=1e3, node.type="red",
                       label.x="Original data-set", label.y="Noise-added data-set")

