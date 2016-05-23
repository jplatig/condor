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
#' @param tick.min.prop minimum proportion of genes/TFs a community must contain to have a visible tick mark, used to prevent overlap of axis labels
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
  reds <- as.character(condor.object$red.memb[order(condor.object$red.memb[,2]),1])
  blues <- as.character(condor.object$blue.memb[order(condor.object$blue.memb[,2]),1])
  if (add.color) {
    # convert edge lists to adjacency matrices (n reds x m blues)
    # plotting colors requires a complete edgelist
    adj <- get.adjacency(condor.object$G, attr="weight", sparse=FALSE)
    # reorder reds according to community membership
    adj <- adj[reds,]
    # reorder blues according to community membership
    adj <- adj[,blues]
    all.links <- melt(adj)
  } else {
    all.links <- cbind(data.frame(get.edgelist(condor.object$G)),get.edge.attribute(condor.object$G,"weight"))
    colnames(all.links) <- c("Var1","Var2","value")
  }
  
  red.sizes <- as.vector(table(condor.object$red.memb[,2]))
  cum.red.sizes <- cumsum(red.sizes)
  blue.sizes <- as.vector(table(condor.object$blue.memb[,2]))
  cum.blue.sizes <- cumsum(blue.sizes)

  ncom <- max(condor.object$blue.memb[,2])
  # reorder by community assignment
  all.links$Var1 <- factor(all.links$Var1,levels=unique(reds))
  all.links$Var2 <- factor(all.links$Var2,levels=unique(blues))
  
  # get axis tick mark locations
  xbreaks <- blues[c(1, cum.blue.sizes[-length(cum.blue.sizes)] + 1)]
  ybreaks <- reds[c(1, cum.red.sizes[-length(cum.red.sizes)] + 1)]

  p <- ggplot() +
    geom_tile(data=all.links, aes(Var2, Var1, fill=value))
  # plot heatmap
  p <- p + xlab(xlab) + 
    ylab(ylab) +
    ggtitle(main) +
    theme(legend.position=legend.position,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  if (ticks=="coms") {
    if (is.null(ylabels)) {
      ylabels <- unique(as.character(sort(condor.object$red.memb[,2])))
      # remove y-axis tick marks that may cause overcrowding
      min.sep <- max(cum.red.sizes) * tick.min.prop
      last.tick.set <- 1
      rm.idx <- c()
      for (i in 1:(length(red.sizes)-1)) {
        if (sum(red.sizes[last.tick.set:i]) > min.sep) {
          last.tick.set <- i + 1
        } else {
          # remove next tick
          rm.idx <- c(rm.idx, i + 1)
        }
      }
      ylabels[rm.idx] <- ""
    }
    if (is.null(xlabels)) {
      xlabels <- unique(as.character(sort(condor.object$blue.memb[,2])))
      # remove x-axis tick marks that may cause overcrowding
      min.sep <- max(cum.blue.sizes) * tick.min.prop
      last.tick.set <- 1
      rm.idx <- c()
      for (i in 1:(length(blue.sizes)-1)) {
        if (sum(blue.sizes[last.tick.set:i]) > min.sep) {
          last.tick.set <- i + 1
        } else {
          # remove next tick
          rm.idx <- c(rm.idx, i + 1)
        }
      }
      xlabels[rm.idx] <- ""
    }
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
  if (transpose) {
    p <- p + coord_flip()
  }
  if (add.color) {
    if (transpose) {
      adj <- t(adj)
    }
    
    combo.grob <- ggplotGrob(p)
    # get subset of adj that contains within-community links
    row.idx <- c(1, cum.red.sizes[-length(cum.red.sizes)] + 1)
    col.idx <- c(1, cum.blue.sizes[-length(cum.blue.sizes)] + 1)
    row.indices <- c()
    col.indices <- c()
    
    # get fill values of  ggplot internals
    rect <- combo.grob$grobs[[4]]$children[[2]]
    fills <- rect$gp$fill
    # reorder matrix
    d <- data.frame(x=as.numeric(rect$x),y=as.numeric(rect$y))
    d$n <- 1:nrow(d)
    d <- d[with(d,order(x,y)),]
    fills <- matrix(fills[d$n],nrow=nrow(adj))
    
    if (is.null(colors)) {
      # cool tutorial http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
      getPalette <- colorRampPalette(brewer.pal(8, "Set1"))
      colors <- getPalette(ncom)
    }
    # convert adjacency matrix values to hex colors
    ConvertMatToHex <- function(mat, xs, ys, color) {
      max.val <- max(mat)
      mat <- mat[xs, ys]
      if (is.null(dim(mat))) {
        sapply(mat,function(x){
          tint <- toupper(format(as.hexmode(round(x/max.val*255)),2))
          paste0(color,tint)
        }) 
      } else {
        apply(mat,2,function(x){
          tint <- toupper(format(as.hexmode(round(x/max.val*255)),2))
          paste0(color,tint)
        })
      }
    }
    for (i in 1:ncom) {
      xs <- row.idx[i]:cum.red.sizes[i]
      ys <- col.idx[i]:cum.blue.sizes[i]
      fills[xs, ys] <- ConvertMatToHex(adj, xs, ys, colors[i])
    }
    # reorder
    fills <- as.vector(fills)[match(1:nrow(d),d$n)]
    combo.grob$grobs[[4]]$children[[2]]$gp$fill <- fills
    plot.new()
    grid.draw(combo.grob)
  } else {
    p
  }
}