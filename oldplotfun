writeOnPlot <- function(txt) {
  plot.new()
  text(0.5, 0.5, txt)
}

scatterplot <- function(data, names, d1, d2, cor=TRUE) {
  plot(data[,d1], data[,d2], xlab=names[d1], ylab=names[d2],
       bty="n", mar=c(0,0,0,0), pch=20)
  if(cor) {
    yrange <- range(data[,d2])
    xrange <- range(data[,d1])
    text(x = 0.1*diff(xrange) + xrange[1], 
         y = 0.9*diff(xrange) + yrange[1],
         labels = paste("r =", signif(cor(data[,d2],data[,d1]),4)))    
  }
}

lineplot <- function(y, xlab=NULL, ylab=NULL, legend.title=NULL) {
  # x - vector; y - vector or matrix
  x <- 1:ncol(y)
  xrange <- range(x)
  yrange <- range(y)
  plot(xrange, yrange, type="n", xlab=xlab, ylab=ylab)
  colors <- rainbow(nrow(y))
  linetype <- c(1:nrow(y))
  plotchar <- linetype + 18
  for (i in 1:nrow(y)) { 
    lines(x, y[i,], type="b", lwd=1.5,
          lty=linetype[i], col=colors[i], pch=plotchar[i]) 
  }
  legend("topleft", legend=rownames(y), cex=0.8, col=colors,
         pch=plotchar, lty=linetype, title=legend.title)
}


## multiplot of ggplots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
