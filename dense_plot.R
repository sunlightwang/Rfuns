source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")

dens_plot <- function(x, y, xlab="x", ylab="y") {
  df <- data.frame(x = x, y = y)
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + scale_fill_gradientn(colours=c("blue","purple", "red"), name = "dens", trans = "log10") + 
  xlab(xlab) + ylab(ylab) + theme_Publication()
}

#require(cowplot)
#  p1 <- dens_plot(log2(cpm[,1]+1), log2(cpm[,2]+1), colnames(cpm)[1], colnames(cpm)[2])
#  p2 <- dens_plot(log2(cpm[,4]+1), log2(cpm[,5]+1), colnames(cpm)[4], colnames(cpm)[5])
#  p3 <- dens_plot(log2(cpm[,1]+1), log2(cpm[,4]+1), colnames(cpm)[1], colnames(cpm)[4])
#  p4 <- dens_plot(log2(cpm[,2]+1), log2(cpm[,5]+1), colnames(cpm)[2], colnames(cpm)[5])
#  plot_grid(p1, p2, p3, p4, ncol=2, labels=c("A", "B", "C", "D"))


# ggplot(df, aes(x=x, y=y, color=class)) + geom_point(size=0.5) + ggtitle("title") + 
#   geom_text_repel(aes(x=x, y=y, color=class, label=names), size=2, segment.size=0.1) + 
#   theme_Publication() 
