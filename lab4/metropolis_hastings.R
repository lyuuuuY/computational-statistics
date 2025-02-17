densityTarget <- function(x) (1 / 14400) * ifelse(x > 0, 120 * x^5 * exp(-x), 0)
densityProposal <- function(x, y) dnorm(x, mean = y, sd = 0.1)
randomProposal <- function(x) rnorm(n = 1, mean = x, sd = 0.1)


#' Metropolis-Hastings Algorithm
#'
#' @param n integer.
#'   The number of samples to generate.
#' @param dtarget function.
#'   The target distribution's density function. 
#'   Must be a function with signature x -> dtarget(x).
#' @param dprop function.
#'   Density function for the proposed conditional distribution. 
#'   Must be a function with signature x, y -> dprop(x, y), where y is the 
#'   conditioning event and x is a value in the support.
#' @param rprop function.
#'   Random generator for the proposed conditional distribution. 
#'   Must be a function with signature x -> rprop(x), where x is the 
#'   conditioning event.
#' @param x0 numeric.
#'   Initial value for the sampling process. Must be in the support of `dprop`.
#'
#' @return numeric vector.
#'   A vector of `n` samples drawn from the target distribution.
rmh <- function(n, densityTarget, densityProposal, randomProposal, x0) {
  
  draws <- numeric(n)
  new_x <- x0
   
  for (i in 1:n) {
    
    last_x <- new_x
    
    # Sample a candidate from the proposal distribution.
    y <- randomProposal(last_x)
    
    # Compute the MH ratio.
    R <- (densityTarget(y) * densityProposal(last_x, y)) / (densityTarget(last_x) * densityProposal(y, last_x))
    
    # Decide new sample
    new_x <- ifelse(
      R >= 1, 
      y,
      ifelse(runif(1) < R, y, last_x)
      )
    
    draws[i] <- new_x
  }
  
  return(draws)
}


x <- rmh(100000, densityTarget, densityProposal, randomProposal, 0.1)
hist(x, probability = TRUE, breaks = 20, col = "skyblue", main = "")
curve(densityTarget(x), add = TRUE, col = "red", lwd = 2)

