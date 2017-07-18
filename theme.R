library(scales)
library(ggplot2)
library(grid)
library(ggthemes)
library(gridExtra)

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.line.x = element_line(color="black"),
           axis.line.y = element_line(color="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.4, "cm"),
           legend.margin = margin(l=0.1, unit="cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
}

theme_graphOnly <- function(base_size=14, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
           text = element_text(),
           axis.line=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = "#f0f0f0"),
           panel.border = element_rect(colour = NA),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.4, "cm"),
           legend.margin = margin(l=0.1, unit="cm"),,
           legend.title = element_text(face="italic")
   ))
}


scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication", manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#a0a013","#dd4444")), ...)
}

scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication", manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#a0a013","#dd4444")), ...)
}

scale_colour_Rainbow <- function(size = 10, ...){
  discrete_scale("colour","Rainbow", manual_pal(values = rainbow(size)), ...)
}
