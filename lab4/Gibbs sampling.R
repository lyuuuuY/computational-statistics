w  <- 1.999
xv <- seq(-1, 1, by=0.01) * 1/sqrt(1-w^2/4)  
plot(xv, xv, type="n", xlab=expression(x[1]), ylab=expression(x[2]), las=1)
# ellipse
lines(xv, -(w/2)*xv-sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=8)
lines(xv, -(w/2)*xv+sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=8)

gibbs_sampling <- function(w,n) {
  x1 <- 0
  x2 <- 0
  samples <- matrix(0, nrow = n, ncol = 2)
  
  for (i in 1:n) {

    x1_min <- -(w/2) * x2 - sqrt(1-(1-w^2/4)*x2^2)
    x1_max <- -(w/2) * x2 + sqrt(1-(1-w^2/4)*x2^2)
    x1 <- runif(1, x1_min, x1_max)
    
    x2_min <- -(w/2) * x1 - sqrt(1-(1-w^2/4)*x1^2)
    x2_max <- -(w/2) * x1 + sqrt(1-(1-w^2/4)*x1^2)
    x2 <- runif(1, x2_min, x2_max)
    
    samples[i, ] <- c(x1, x2)
  }
  
  return(samples)
}

set.seed(12345)
samples <- gibbs_sampling(1.999,1000)
mean(samples[,1] > 0)
points(samples[,1], samples[,2], pch = 20, cex = 0.5, col = "blue")

 
gibbs_sampling_U <- function(w, n) {
  U1 <- 0
  U2 <- 0
  samples_U <- matrix(0, nrow = n, ncol = 2)
  
  for (i in 1:n) {
    
    molecular1 <- 4 - (2 + w)*U2^2
    rangeU1 <- sqrt(molecular1 / (2 - w))
    U1_min <- -rangeU1
    U1_max <-  rangeU1
    U1 <- runif(1, U1_min, U1_max)
    
    molecular2 <- 4 - (2 - w)*U1^2
    rangeU2 <- sqrt(molecular2 / (2 + w))
    U2_min <- -rangeU2
    U2_max <-  rangeU2
    U2 <- runif(1, U2_min, U2_max)
    
    samples_U[i, ] <- c(U1, U2)
  }
  

  return(samples_U)
}
samples_U <- gibbs_sampling_U(1.999,1000)
samples_X <- cbind( (samples_U[,1] + samples_U[,2]) / 2,
                    (samples_U[,2] - samples_U[,1]) / 2 )


u1v <- seq(-1, 1, by = 0.01) * sqrt(4 / (2 - w))
u2_plus  <-  sqrt( pmax(0, 4 - (2 - w)*u1v^2 ) / (2 + w) )
u2_minus <- -u2_plus
plot(u1v, u2_plus, type="n", xlab=expression(u[1]), ylab=expression(u[2]), 
     asp = 1, las = 1, 
     main = bquote("Ellipse in (u1,u2):"~(2-w)*u[1]^2 + (2+w)*u[2]^2==4))
lines(u1v,  u2_plus,  lwd=2, col=8)
lines(u1v,  u2_minus, lwd=2, col=8)

points(samples_U[,1], samples_U[,2], pch = 20, cex = 0.5, col = "red")
