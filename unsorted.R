# test bimodality
library(diptest)
dip.test(c(rnorm(1000), rnorm(1000,3)))
