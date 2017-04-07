library(matrixStats)

G1 <- rnorm(1000)
G2 <- rnorm(1000)*20
GP <- rbind(G1, G2)
plot(t(GP))
center.scale <- function(data)  (data - rowMeans(data)) / sqrt(rowVars(data))
GP.scaled <- center.scale( GP )
plot(t(GP.scaled))
pca <- prcomp(t(GP.scaled))
# pca$x = t(GP.scaled) %*% pca$rotation

pca.rotated <- t(GP.scaled) %*% pca$rotation
pca.rotated.scaled <- t( t(pca.rotated) / pca$sdev )
plot(pca.rotated.scaled)
abline(v=0,h=0,lty=2,col="grey80")

cart2pol <- function(data) {
  do.call(rbind, sapply(1:nrow(data), function(x) {
    r <- sqrt(data[x,1]^2 + data[x,2]^2)
    theta <- atan(data[x,2] / data[x,1])
    if(data[x,2]<0 & data[x,1]<0) {theta=theta-pi}
    if(data[x,2]>0 & data[x,1]<0) {theta=theta+pi}
    return(cbind(r, theta))
  }, simplify=F))
}

GP.pole <- cart2pol(pca.rotated.scaled)
plot(GP.pole)
hist(GP.pole[,1])
