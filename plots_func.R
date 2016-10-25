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
  r <- abs(cor(x[!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)], y[!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)]))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.68/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (r+1)/2 )
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

