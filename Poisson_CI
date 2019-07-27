
Poisson_CI <- function(N1=10000, N2=10000, CI=0.95, plot=T) {
  # N1 <- 10000
  # N2 <- 10000
  ci_p <- qnorm((1-CI)/2, lower.tail = F)
  x <- 1:min(N1,N2)
  y <- 1:min(N1,N2)
  p1 <- x/N1
  p2 <- y/N2
  delta_x <- sqrt((1 - p1) / x) * log10(exp(1))
  delta_y <- sqrt((1 - p2) / y) * log10(exp(1))
  delta_M <- 2 * sqrt(delta_x ^ 2 * delta_y ^ 2 / (delta_x ^ 2 + delta_y ^ 2) ) 
  mu <- (x + y) / sqrt(2)
  x_p1 <- log10(mu) / sqrt(2) - (ci_p * delta_M - log10(N1/N2)) / sqrt(2) 
  y_p1 <- log10(mu) / sqrt(2) + (ci_p * delta_M - log10(N1/N2)) / sqrt(2) 
  x_p2 <- log10(mu) / sqrt(2) - (-ci_p * delta_M - log10(N1/N2)) / sqrt(2) 
  y_p2 <- log10(mu) / sqrt(2) + (-ci_p * delta_M - log10(N1/N2)) / sqrt(2) 
  if(plot) {
    plot(x_p1, y_p1, type="l", lty=2, xlim=c(-1,log10(N1)-0.5), ylim=c(-1,log10(N2)-0.5))
    lines(x_p2, y_p2, type="l", lty=2)
  } else {
    return(rbind(data.frame(x=x_p1, y=y_p1, name="up"), data.frame(x=x_p2, y=y_p2, name="down")))
  }
}
