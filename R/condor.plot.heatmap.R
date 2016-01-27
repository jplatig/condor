#' Plot weighted adjacency matrix with links grouped by community
#' 
#' This function will generate the network link 'heatmap' for a weighted network
#' @param condor.object output of either \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @param ticks choose to display row and column names, red and blue community numbers or no axis ticks
#' @param add.color if TRUE, fill in communities with solid colors
#' @param colors character vector of colors of length equal to the number of communities
#' in the network
#' @param log if TRUE, transform the color gradient to log-scale.
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main a title for the plot
#' @param tick.min.prop minimum proportion of genes/TFs a community must contain to have a visible tick mark, used to prevent overlap
#' @param xlabels optional character vector of length equal to the number of communities to manually mark the x-axis community ticks
#' @param ylabels optional character vector of length equal to the number of communities to manually mark the y-axis community ticks
#' @param legend.position side the legend is placed
#' @param transpose if TRUE, plot reds in columns rather than rows
#' @return produces a ggplot object if \code{color} is FALSE, else no return value
#' @examples
#' data(small1976)
#' condor.object <- create.condor.object(small1976)
#' condor.object <- condor.cluster(condor.object, project=FALSE)
#' condor.plot.heatmap(condor.object)
#' @importFrom grid addGrob getGrob grid.draw
#' @import ggplot2
#' @import RColorBrewer
#' @export
#'  
condor.plot.heatmap = function(condor.object, ticks=c("names", "coms", "none"), add.color=FALSE,
                               colors=NULL, log=FALSE, xlab="blues", ylab="reds",
                               main="", tick.min.prop=0.05, xlabels=NULL, ylabels=NULL,
                               legend.position="right", transpose=FALSE){
  ticks <- match.arg(ticks)
  # set edges in unweighted graph to have weight of 1
  if (is.null(get.edge.attribute(condor.object$G, "weight"))) {
    condor.object$G <- set.edge.attribute(condor.object$G, "weight", value=1)
  }
  # convert edge lists to adjacency matrices (n reds x m blues)
  adj <- get.adjacency(condor.object$G, attr="weight", sparse=FALSE)
  if (transpose) {
    adj <- t(adj)
    tmp <- condor.object$blue.memb
    condor.object$blue.memb <- condor.object$red.memb
    condor.object$red.memb <- tmp
    if (is.null(xlab)) {
      xlab <- "reds"
    }
    if (is.null(ylab)) {
      ylab <- "blues"      
    }
  }
  # reorder reds according to community membership
  reds <- as.character(condor.object$red.memb[order(condor.object$red.memb[,2]),1])
  adj <- adj[reds,]
  # reorder blues according to community membership
  blues <- as.character(condor.object$blue.memb[order(condor.object$blue.memb[,2]),1])
  adj <- adj[,blues]
  rowsep <- cumsum(as.vector(table(condor.object$red.memb[,2])))
  colsep <- cumsum(as.vector(table(condor.object$blue.memb[,2])))
  ncom <- max(condor.object$blue.memb[,2])
  all.links <- melt(adj)
  
  # get axis tick mark locations
  xbreaks <- colnames(adj)[c(1, colsep[-length(colsep)] + 1)]
  ybreaks <- rownames(adj)[c(1, rowsep[-length(rowsep)] + 1)]
  if (!transpose) {
    if (is.null(xlabels)) {
      xlabels <- unique(as.character(sort(condor.object$blue.memb[,2])))
      # remove x-axis tick marks that may cause overcrowding
      idx <- which(diff(colsep) < max(colsep) * tick.min.prop)[1]
      if (!is.na(idx)) {
        xlabels[seq(idx + 2, length(xbreaks), 2)] <- ""  
      }
    }
    if (is.null(ylabels)) {
      ylabels <- unique(as.character(sort(condor.object$blue.memb[,2])))
      # remove y-axis tick marks that may cause overcrowding
      row.diff <- diff(rowsep)
      row.rle <- rle(row.diff)
      idx <- which(row.rle$values < max(rowsep) * tick.min.prop)
      if (length(idx) > 0) {
        rm.idx <- c()
        for(i in idx) {
          csum <- sum(row.rle$len[1:(idx-1)])
          rm.idx <- c(rm.idx, seq(csum + 2, csum + row.rle$lengths[i] + 1, 2))
        }
        ylabels[rm.idx] <- ""  
      }
    }
  } else {
    if (is.null(ylabels)) {
      ylabels <- unique(as.character(sort(condor.object$blue.memb[,2])))
      # remove x-axis tick marks that may cause overcrowding
      idx <- which(diff(rowsep) < max(rowsep) * tick.min.prop)[1]
      if (!is.na(idx)) {
        ylabels[seq(idx + 2, length(ybreaks), 2)] <- ""  
      }
    }
    if (is.null(xlabels)) {
      xlabels <- unique(as.character(sort(condor.object$blue.memb[,2])))
      # remove y-axis tick marks that may cause overcrowding
      row.diff <- diff(colsep)
      row.rle <- rle(row.diff)
      idx <- which(row.rle$values < max(colsep) * tick.min.prop)
      if (length(idx) > 0) {
        rm.idx <- c()
        for(i in idx) {
          csum <- sum(row.rle$len[1:(idx-1)])
          rm.idx <- c(rm.idx, seq(csum + 2, csum + row.rle$lengths[i] + 1, 2))
        }
        xlabels[rm.idx] <- ""  
      }
    }
  }
  
  p <- ggplot() +
    geom_tile(data=all.links, aes(Var2, Var1, fill=value))
  # plot heatmap
  p <- p + xlab(xlab) + 
    ylab(ylab) +
    ggtitle(main) +
    theme(legend.position=legend.position)
  if (ticks=="coms") {
    p <- p + scale_y_discrete(breaks=ybreaks, labels=ylabels) +
      scale_x_discrete(breaks=xbreaks, labels=xlabels)
  }
  if (ticks=="names") {
    p <- p + theme(axis.text.x=element_text(angle = 45, hjust = 1))
  }
  if (ticks=="none") {
    p <- p + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  }
  if (log) {
    p <- p + scale_fill_gradient(name="log10(edge weight)", low="white", high="gray20", trans="log10", na.value="white")
  } else {
    p <- p + scale_fill_gradient(name="Edge weight", low="white", high="gray20", na.value="white")
  }
  if (add.color) {
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
      p2 <- ggplot() +
        geom_tile(data=com.links, aes(Var2, Var1, fill=value))
      if (ticks=="coms") {
        p2 <- p2 + scale_y_discrete(breaks=ybreaks, labels=ylabels) +
          scale_x_discrete(breaks=xbreaks, labels=xlabels)
      }
      if (ticks=="names") {
        p2 <- p2 + theme(axis.text.x=element_text(angle = 45, hjust = 1))
      }
      if (ticks=="none") {
        p2 <- p2 + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
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
    plot.new()
    grid.draw(combo.grob)
  } else {
    p
  }
}