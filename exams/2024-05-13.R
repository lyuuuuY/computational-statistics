
#  --- ASSIGNMENT 1 ---

# FUNction, Gradient and Hessian.

FUN <- function(x, y) -x^2 - x^2*y^2 - 2*x*y + 2*x + 2

G <- function(X) {
  x <- X[1]
  y <- X[2]
  
  matrix(
    c(-2*x - 2*x*y^2 - 2*y + 2, -2*x^2*y - 2*x),
    nrow = 2,
    ncol = 1
  )
}

H <- function(X) {
  x <- X[1]
  y <- X[2]
  
  matrix(
    c(-2 - 2*y^2, -4*x*y - 2, -4*x*y - 2, -2*x^2),
    nrow = 2,
    ncol = 2
  )
}


# >>> Question 1.1

#' Newton optimizer.
#'
#' @param x0 vector. Initial guess
#' @param FUN function. Function to be optimized
#' @param G function. Gradient
#' @param H function. Hessian
#' @param tol float. Stopping tolerance
#'
#' @returns vector. Optimizer.
mynewton <- function(x0, FUN, G, H, tol = 0.0001) {
  
  # IDEA: at each iteration, approximate the real optimum using the 
  # second order Taylor approximation.
  
  delta <- 1e-8
  
  # Step for Newton Method.
  step <- function(x) x - solve(H(x)) %*% G(x)
  
  # Perform initial step to initialize variables.
  x <- x0
  xnew <- step(x)
  conv <- norm((xnew - x) / (x + delta))
  
  while (conv > tol) {
    x <- xnew
    xnew <- step(x)
    conv <- norm((xnew - x) / (x + delta))
    
  }
  return (xnew)
}

# Contour plot for sanity check.
x <- y <- seq(-3, 3, length.out = 100)
z <- outer(x, y, FUN)


# >>> Question 1.2
x0_1 <- c(2, 0)
x0_2 <- c(-1, -2)
x0_3 <- c(0, 1)
# Using `x0_1` returns (1, -1) while `x0_2` and `x0_3` return (0, 1) which
# suggest both (x, y) correspond to points where the gradient vanishes due
# to either saddle points or maxima.


# >>> Question 1.3
# To extend the algorithm for finding a local maximum, check the Hessian matrix 
# H(x) at each iteration and change the Newton direction accordingly.



# --- ASSIGNMENT 2 ---

# We are interested to generate draws of a random variable X with this density 
# using rejection sampling.

# >>> Question 2.1
# Good Envelope: Exponential distribution with lambda = 1. 
# This distribution is also only defined for values where the target density
# exists and will effectively majorize the target density for all values of x
# while not wasting too much space that can result in poor accepting rate.

# Bad Envelope: A normal distribution.
# Normal distributions are define for the entire real line, thus causing
# unnecessary sample values since the target is only defined for positive.

# Another Bad Envelope: Uniform distribution.
# Uniform distribution is only suitable when the target support is finite.


# >>> Question 2.2

# Define densities.
# The proposed variable Y is chosen so that we
# can easily generate realizations of it and so that its density g_Y can be scaled
# to majorize p_X using some constant c; i.e., c * g_Y > p_X

dtarget <- function(x) ifelse(x >= 0, (5 / 6) * exp(-x) * (1 + cos(2*x)), 0)
denv_good <- function(x) 2 * dexp(x)
denv_bad <- function(x) dunif(x)

# Define random number generator for proposed envelope distribution.
renv <- function() rexp(1)

# Evaluate densities in [0, 5]
x <- seq(0, 5, length.out = 50)
y_target <- dtarget(x)
y_good <- denv_good(x)
y_bad <- denv_bad(x)

# Plot densities
plot(x, y_target, type="l", col = "black", lwd = 1)
lines(x, y_good, col = "blue")
lines(x, y_bad, col = "orange")

legend("topright", legend = c("Target Density", "Good Envelope", "Bad Envelope"), 
       col = c("black", "blue", "orange"), lwd = 2)



# >>> Question 2.3

#' Factory for rejection sampler
#'
#' @param renv. Function. Random number generator for proposed/envelope distribution
#' @param denv. Function. Density for proposed/envelope distribution
#' @param dtarget. Function. Density for target distribution
#'
#' @returns Function. A random number generator
make_rejection_sampler <- function(renv, denv, dtarget) {
  
  rejection_sampler <- function(n) {
    
    sample <- numeric(n)
    count <- 0
    it <- 0
    
    while (count < n) {
      
      # Generate y from the distribution with density function gY.
      y <- renv()
      
      # Generate u from a uniform (0,1) distribution.
      u <- runif(1)
      
      # Sampling ratio
      ratio <- dtarget(y) / denv(y)
      
      if (ratio >= u) {
        # Accepted sample.
        sample[count] <- y
        count <- count + 1
      }
      
      it <- it + 1
      
    }
    
    print(n / it)
    return(sample)
  }
  return(rejection_sampler)
}

n <- 10000
rej_sampler <- make_rejection_sampler(renv, denv_good, dtarget)
sample <- rej_sampler(n)
hist(sample, breaks = 30, probability = TRUE)
curve(dtarget(x), add = TRUE, col = "red", lwd = 2)


# >>> Question 2.4

# Based on the sample generated in Part 2.3, determine estimates for P (X > 1) and E(X).
# Determine some measure of uncertainty for these estimates

p_x_greater_than_1 <- sum(sample > 1) / n
expected_value <- mean(sample)