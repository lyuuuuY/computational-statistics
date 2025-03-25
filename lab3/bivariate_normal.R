# Generate a two-dimensional random vector [X, Y] that has a two-dimensional normal distribution
# with mean vector [0, 0] and and covariance matrix [[0.6, 0], [0, 0.6]] using 
# the Box-Muller method and using only runif as random number generator.

rnomalbm <- function(n, mean = 0, sd = 1) {
  
  U <- runif(n)
  V <- runif(n, max = 2*pi)
  
  X <- sqrt(-2 * log(U)) * cos(V)
  Y <- sqrt(-2 * log(U)) * sin(V)
  
  Z <- cbind(X, Y)
  mean + sd * Z
}
