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