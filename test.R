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

  ####
pol2cart <- function(data) {
  do.call(rbind, sapply(1:nrow(data), function(x) {
    xx <- data[x,1] * cos(data[x,2])
    yy <- data[x,1] * sin(data[x,2])
    return(cbind(xx, yy))
  }, simplify=F))
}

r <- rnorm(1000,1,0.2)
th <- runif(1000, -pi, pi)
circ <- t(pol2cart(cbind(r,th)))
plot(t(circ))
circ.scaled <- center.scale( circ )
plot(t(circ.scaled))
pca <- prcomp(t(circ.scaled))
pca.rotated <- t(circ.scaled) %*% pca$rotation
pca.rotated.scaled <- t( t(pca.rotated) / pca$sdev )
plot(pca.rotated.scaled)
abline(v=0,h=0,lty=2,col="grey80")
GP.pole <- cart2pol(pca.rotated.scaled)
plot(GP.pole)
hist(GP.pole[,1])
