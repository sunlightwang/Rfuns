panel.hist <- function(x, breaks=20, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, breaks=breaks, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x[!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)], y[!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)])
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.68/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (abs(r)+1)/2 )
}

panel.smooth <- function(x, y, col = par("col"), pch = par("pch"), cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  smoothScatter(x, y, pch = pch, colramp = colorRampPalette(c("white","lightgreen","blue","purple"), bias=1.2), cex = cex,
                nrpoints = 100, col=col, nbin = 256, add=TRUE, ...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), col = col.smooth, ...)
  abline(0, 1, lty=2, lwd=0.1, col=col.smooth)
}


panel.diagonal <- function(x, y, col = par("col"), pch = par("pch"), cex = par("cex"), col.diagonal = "blue", ...) 
{
  points(x, y, pch = pch, cex = cex, col=col,  ...)
  abline(v=0, h=0, lty=2, lwd=0.1, col=col.diagonal)
  abline(-1, 1, lty=2, lwd=0.1, col=col.diagonal)
  abline(0, 1, lty=2, lwd=0.1, col=col.diagonal)
  abline(1, 1, lty=2, lwd=0.1, col=col.diagonal)
}

#### old functions
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
