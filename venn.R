draw.venn <- function(set1, set2, group_names=NULL, doWeights = TRUE) {
  library(Vennerable)
  data <- list(set1, set2)
  names(data) <- group_names
  V <- Venn(data)
  plot(V, doWeights = doWeights)
}
draw.venn3 <- function(set1, set2, set3,  group_names=NULL, doWeights = TRUE) {
  library(Vennerable)
  data <- list(set1, set2, set3)
  names(data) <- group_names
  V <- Venn(data)
  plot(V, doWeights = doWeights)
}
draw.venn4 <- function(set1, set2, set3, set4,  group_names=NULL, doWeights = FALSE) {
  library(Vennerable)
  data <- list(set1, set2, set3, set4)
  names(data) <- group_names
  V <- Venn(data)
  plot(V, doWeights = doWeights, type="ellipses")
}
