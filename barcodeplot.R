library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

barcodeplot <- function(data, sample1="rep1", sample2="rep2", colorby=NULL) {
  ## selecte data
  colorby <- match.arg(colorby, c("size", colnames(data)))
  data <- data.frame(data, size=factor(nchar(rownames(data)), levels=c(1,3,5,7,9)))
  data0 <- data.frame(x=data[,sample1], y=data[,sample2], col=data[,colorby], row.names=rownames(data))
  data0 <- data0[apply(data0[,c("x","y")], 1, sum)>0, ]
  data1 <- data0
  data1$x[data1$x==0] <- 0.1
  data1$y[data1$y==0] <- 0.1
  #data2 <- data0[apply(data0[,c("x","y")]>0, 1, all), c("x","y")]
  data2 <- data0[apply(data0[,c("x","y")]>0, 1, any), c("x","y")]
  ## both non-zero correlation 
  r <- signif(cor(data2, method="spearman")[2,1],3)
  ## one-zero range
  r1 <- range(data0[data0[,"y"]==0, "x"])
  r2 <- range(data0[data0[,"x"]==0, "y"])
  ## polygon
  datapoly.x <- data.frame(x=c(r1[1]*0.8, r1[2]*1.25, r1[2]*1.25,  r1[1]*0.8), y=c(0.08,0.08,0.125,0.125))
  datapoly.y <- data.frame(y=c(r2[1]*0.8, r2[2]*1.25, r2[2]*1.25,  r2[1]*0.8), x=c(0.08,0.08,0.125,0.125))
  ## one-zero percentage
  p.x <- signif(sum(data0[,"y"]==0)/nrow(data0) * 100, 3)
  p.y <- signif(sum(data0[,"x"]==0)/nrow(data0) * 100, 3)
  data.text <- data.frame(x=c(0.1, r1[2]*1.3), y=c(r2[2]*1.3, 0.1), text=paste0(c(p.y,p.x),"%"))
  data.text <- rbind(data.text, data.frame(x=1,y=max(data1$y)*0.8,text=paste0("Rs=",r)))
  ## plots
  if(!is.null(colorby)) {
    p <- ggplot(data1) + geom_point(aes(x,y,color=col)) + labs(color=colorby) 
  } else {
    p <- ggplot(data1) + geom_point(aes(x,y))
  }
  p <- p + scale_x_continuous(limits=c(0.08, max(data1$x)*1.25), trans='log10',breaks=c(0.1, 1, 10,100,1000,10000), labels=c(0,1,10,100,1000,10000))  +  
    scale_y_continuous(limits=c(0.08, max(data1$y)*1.25), trans='log10', breaks=c(0.1, 1, 10,100,1000,10000), labels=c(0,1,10,100,1000,10000)) + 
    xlab(sample1) + ylab(sample2)
  p <- p + geom_polygon(data=datapoly.x, mapping=aes(x,y), fill="grey50",alpha=0.3) + 
    geom_polygon(data=datapoly.y, mapping=aes(x,y), fill="grey50",alpha=0.3)
  p + geom_text(aes(x,y,label = text), data=data.text, hjust=0, vjust=0)
}
