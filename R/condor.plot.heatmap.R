#' Plot weighted adjacency matrix with links grouped by community
#' 
#' This function will generate the network link 'heatmap' for a weighted network
#' @param condor.object output of either \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @param community.ticks if TRUE, plot community numbers in x and y axes rather than column and row names
#' @param color if TRUE, fill in communities with solid colors
#' @param colors character vector of colors of length equal to the number of communities
#' in the network
#' @param log if TRUE, transform the color gradient to log-scale.
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @return produces a ggplot object if color is FALSE, else no return value
#' @examples
#' data(small1976)
#' condor.object <- create.condor.object(small1976)
#' condor.object <- condor.cluster(condor.object, project=FALSE)
#' condor.plot.heatmap(condor.object)
#' @import ggplot2
#' @import grid
#' @import RColorBrewer
#' @export
#'  
condor.plot.heatmap = function(condor.object, community.ticks=FALSE, color=FALSE, colors=NULL, log=FALSE, xlab="blues", ylab="reds"){
  # set edges in unweighted graph to have weight of 1
  if (is.null(get.edge.attribute(condor.object$G, "weight"))) {
    condor.object$G <- set.edge.attribute(condor.object$G, "weight", value=1)
  }
  # convert edge lists to adjacency matrices (n reds x m blues)
  adj <- get.adjacency(condor.object$G, attr="weight", sparse=FALSE)
  # reorder reds according to community membership
  reds <- as.character(condor.object$red.memb[order(condor.object$red.memb[,2]),1])
  adj <- adj[reds,]
  # reorder blues according to community membership
  blues <- as.character(condor.object$blue.memb[order(condor.object$blue.memb[,2]),1])
  adj <- adj[,blues]
  rowsep <- cumsum(as.vector(table(condor.object$red.memb[,2])))
  colsep <- cumsum(as.vector(table(condor.object$blue.memb[,2])))
  lab <- unique(as.character(sort(condor.object$blue.memb[,2])))
  ncom <- length(lab)
  all.links <- melt(adj)
  
  p <- ggplot() +
    geom_tile(data=all.links, aes(Var2, Var1, fill=value)) +
    xlab(xlab) + 
    ylab(ylab)
  if (community.ticks) {
    p <- p + scale_y_discrete(breaks=rownames(adj)[c(1, rowsep[-length(rowsep)])], labels=lab) +
      scale_x_discrete(breaks=colnames(adj)[c(1, colsep[-length(colsep)])], labels=lab)
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  if (log) {
    p <- p + scale_fill_gradient(name="log10(edge weight)", low="white", high="gray20", trans="log10", na.value="white")
  } else {
    p <- p + scale_fill_gradient(name="Edge weight", low="white", high="gray20", na.value="white")
  }
  if (color) {
    combo.grob <- ggplotGrob(p)
    # get subset of adj that contains within-community links
    row.idx <- c(1, rowsep[-length(rowsep)] + 1)
    col.idx <- c(1, colsep[-length(colsep)] + 1)
    row.indices <- c()
    col.indices <- c()
    
    if (is.null(colors)) {
      # cool tutorial http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
      getPalette <- colorRampPalette(brewer.pal(8, "Set1"))
      colors <- getPalette(ncom)
    }
    for (i in 1:ncom) {
      m <- matrix(0, nrow=nrow(adj), ncol=ncol(adj))
      rownames(m) <- rownames(adj)
      colnames(m) <- colnames(adj)
      m[row.idx[i]:rowsep[i], col.idx[i]:colsep[i]] <- adj[row.idx[i]:rowsep[i], col.idx[i]:colsep[i]]
      # hack to ensure consistent color scaling
      if (i != ncom) {
        m[nrow(m), ncol(m)] <- max(adj)  
      } else {
        m[1, 1] <- max(adj)
      }
      com.links <- melt(m)
      p2 <- ggplot() + geom_tile(data=com.links, aes(Var2, Var1, fill=value))
      if (community.ticks) {
        p2 <- p2 + scale_y_discrete(breaks=rownames(adj)[c(1, rowsep[-length(rowsep)])], labels=lab) +
          scale_x_discrete(breaks=colnames(adj)[c(1, colsep[-length(colsep)])], labels=lab)
      } else {
        p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
      if (log) {
        p2 <- p2 + scale_fill_gradient(name="log10(edge weight)", low="white", high=colors[i], trans="log10", na.value="white")
      } else {
        p2 <- p2 + scale_fill_gradient(name="edge weight", low="white", high=colors[i], na.value="white")
      }
      p2g <- ggplotGrob(p2)
      # change white to transparent
      p2g$grobs[[4]]$children[[2]]$gp$fill <- gsub("#FFFFFFFF", "#00000000", p2g$grobs[[4]]$children[[2]]$gp$fill)
      if (i == ncom) {
        p2g$grobs[[4]]$children[[2]]$gp$fill[1] <- "#00000000"
      } else {
        p2g$grobs[[4]]$children[[2]]$gp$fill[length(p2g$grobs[[4]]$children[[2]]$gp$fill)] <- "#00000000"
      }
      combo.grob$grobs[[4]] <- addGrob(combo.grob$grobs[[4]],
                                       getGrob(p2g$grobs[[4]],
                                               "geom_rect.rect",
                                               grep=TRUE))
    }
    grid.draw(combo.grob)
  } else {
    p
  }
  
}
