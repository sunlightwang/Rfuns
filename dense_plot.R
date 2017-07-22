dens_plot <- function(x, y, xlab="x", ylab="y") {
  df <- data.frame(x = x, y = y)
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + scale_fill_gradientn(colours=c("blue","purple", "red"), name = "dens", trans = "log10") + 
  xlab(xlab) + ylab(ylab) + theme_Publication()
}
