# Draw the boundaries of the region where X has a uniform distribution.
w  <- 1.8
xv <- seq(-1, 1, by=0.01) * 1/sqrt(1-w^2/4)  # a range of x1-values, where the term below the root is non-negative (compare Lecture 4)
plot(xv, xv, type="n", xlab=expression(x[1]), ylab=expression(x[2]), las=1)
# ellipse
lines(xv, -(w/2)*xv-sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=8)
lines(xv, -(w/2)*xv+sqrt(1-(1-w^2/4)*xv^2), lwd=2, col=8)


# What is the conditional distribution of X1 given X2 and that of X2 given X1?
# Both are uniformly distributed over their corresponding ellipse cross section.


# Write your own code for Gibbs sampling the distribution.
gibbs_sampler <- function(n, w) {
  x1 <- 0
  x2 <- 0
  sample <- matrix(nrow = n, ncol = 2)
  
  for (i in 1:n) {
    x1_min <- -(w/2) * x2 - sqrt(1-(1-w^2/4)*x2^2)
    x1_max <- -(w/2) * x2 + sqrt(1-(1-w^2/4)*x2^2)
    x1 <- runif(1, x1_min, x1_max)
    
    x2_min <- -(w/2) * x1 - sqrt(1-(1-w^2/4)*x1^2)
    x2_max <- -(w/2) * x1 + sqrt(1-(1-w^2/4)*x1^2)
    x2 <- runif(1, x2_min, x2_max)
    
    sample[i,] <- c(x1, x2)
  
  }
  return(sample)
}


S <- gibbs_sampler(1000, w)
points(S[, 1], S[, 2], pch = 16, cex = .2)

