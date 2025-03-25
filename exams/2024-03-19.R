# --- ASSIGNMENT 1 ---

# >>> Question 1.1
# Plot the function f.
x <- seq(-20, 20, length.out = 5000)
f <- function(x) exp(-x^2 / 100) * (2 - 0.1 * cos(3* pi * (x + 0.1)))
plot(x, f(x), type = "l")


# >>> Question 1.2
mysimann <- function(FUN, x0, t0, tol = 1e-6) {
  
  n_iterations <- 1000
  
  # Proposal distribution numbers generator.
  # Here, we simulate a random walk using a uniform U(x-a, x+a) with a = 5.
  rproposed <- function(x) runif(1, x - 1, x + 1)
  
  # Cooling function.
  # Straight line with negative slope.
  alpha <- function(t) -0.5 * t
  
  # Initialize variables.
  # Set `conv` high so the while loop is triggered.
  t <- t0
  xnew <- x0
  conv <- 10000

  while (conv > tol) {
    
    for (i in 1:n_iterations) {
      x <- xnew
      
      # Sample candidate from proposal distribution,
      # i.e., a neighboring solution.
      y <- rproposed(x)
      
      # Delta Energy.
      C <- FUN(y) - FUN(x)
      
      if (C > 0) {
        # Accept candidate (new state)
        xnew <- y
      } else {
        # Compute h.
        h <- exp(FUN(y) - FUN(x) / t)
        
        # Reject/Accept candidate.
        xnew <- ifelse(runif(1) < h, y, x)
      }
    }
    conv <- abs(xnew - x)
    t <- alpha(t)
  }
  
  return(xnew)
  
}


# --- ASSIGNMENT 2 ---

# The task is here to generate a random sample following this bivariate normal distribution.

# >>> Question 2.1
# Implement a Gibbs sampler to generate a random sample of size n for X with this bivariate
# normal distribution.

gibbs <- function(n, x0, ro) {
  
  d <- length(x0)  # d = 2 (bivariate normal)
  sample <- matrix(nrow = n, ncol = d)
  i <- 0
  x <- x0
  
  for (i in 1:n) {
    
    x1 <- x[1]
    x2 <- x[2]
    
    # Sample X1 | X2
    x1new <- rnorm(1, mean = ro * x2, sd = sqrt(1 - ro^2))
    
    # Sample X2 | X1 (using the updated X1)
    x2new <- rnorm(1, mean = ro * x1new, sd = sqrt(1 - ro^2))
    x <- c(x1new, x2new)
    
    sample[i,] <- x
  
  }
  
  return(sample)
}


# >>> Question 2.2
n <- 10000      # Number of samples
x0 <- c(0, 0)  # Initial value

# ro = 0
ro <- 0        # Correlation coefficient
sample1 <- gibbs(n, x0, ro)
plot(
  sample1[, 1], 
  sample1[, 2], 
  xlab = "X1", 
  ylab = "X2", 
  main = "Gibbs Sampler for Bivariate Normal"
)

# ro = 0.998
ro <- 0.998        # Correlation coefficient
sample2 <- gibbs(n, x0, ro)
plot(
  sample2[, 1], 
  sample2[, 2], 
  xlab = "X1", 
  ylab = "X2", 
  main = "Gibbs Sampler for Bivariate Normal"
)



# >>> Question 2.3
# Based on your generated samples, provide an approximation for P (X1 + X2 > 0) 
# in both cases.

sum1 <- sample1[, 1] + sample1[, 2]
prob1 <- sum(sum1 > 0) / n  # Sample fraction estimation.

sum2 <- sample2[, 1] + sample2[, 2]
prob2 <- sum(sum2 > 0) / n  # Sample fraction estimation.


# >>> Question 2.4
#Instead of using Gibbs sampling, suggest an alternative algorithm which can be used to generate
# such a bivariate normal distribution and describe it

# Suggestion: Rejection sampling.
# Description: Rejection sampling uses consists of choosing a random variable Y with density g_Y
# that resembles the target density and from which it is easy to sample from.
# The density `g_y` is known was the envelope density and must satisfy c*g_Y > p_X
# for all x. That is, it acts as an envelope for the target density. In fact,
# how good of an envelope g_Y is determines the quality of the sampling.
# A random sample from Y is accepted only when p_X(y) / (c*g_Y(y)) > U.
# Cons: Rejection sampling works poorly in higher dimensions since even a thin 
# shell between both densities contains most of the volume, so the rejection
# proportion would be high.



