#  ---------------- Factory functions ----------------

#' Creates an alpha scheduler function for step size adjustment.
#'
#' @param factor Numeric. Reduction factor for alpha when no improvement 
#' is observed. Default is 0.1.
#'
#' @return A function that adjusts alpha based on function evaluations.
make_alpha_scheduler <- function(factor = 0.1) {
  
  scheduler <- function(alpha, val, newval) {
       
      if (val < newval) {
        # If val < newval (improvement), do nothing.
        return(alpha)
      } else {
        # Otherwise, reduce ´alpha´.
        return(alpha * factor)
      }
    }
  return(scheduler)
} 

#' Creates a log-likelihood function for logistic regression.
#'
#' @param x: A numeric vector of predictor values.
#' @param y: A numeric vector of binary response values (0 or 1).  
#'
#' @return: A function that computes the log-likelihood given a parameter vector TH.
make_llik <- function(x, y) {
  
  llik <- function(TH) {
    theta_0 <- TH[1]
    theta_1 <- TH[2]
    p <- 1 / (1 + exp(- theta_0 - theta_1 * x))
    log_likelihood <- sum(y * log(p) + (1 - y) * log(1 - p))
    
    return(log_likelihood)
  }
  
  return(llik)
}


#' Creates the gradient function of the log-likelihood function for logistic regression.
#'
#' @param x A numeric vector of predictor values.
#' @param y A numeric vector of binary response values (0 or 1).
#'
#' @return A function that computes the gradient given a parameter vector TH.
make_grad <- function(x, y) {
  
  grad <- function(TH) {
    theta_0 <- TH[1]
    theta_1 <- TH[2]
    
    p <- 1 / (1 + exp(-(theta_0 + theta_1 * x)))
    gradient <- c(sum(y - p), sum((y - p) * x))
    
    return(gradient)
  }
  
  return(grad)
}


#  ---------------- Steepest ascent ----------------


#' Performs steepest ascent optimization with optional step-size adaptation.
#'
#' @param f Function. The objective function to maximize.
#' @param g Function. The gradient function.
#' @param x0 Numeric. Initial point.
#' @param eps Numeric. Convergence threshold. Default is 1e-8.
#' @param alpha.0 Numeric. Initial step size. Default is 1.
#' @param alpha.scheduler Function. Step-size adjustment strategy. 
#' Default is NULL (constant step size).
#'
#' @return A data frame tracking the iteration history.
steepestasc <- function(
    f, 
    g, 
    x0, 
    eps = 1e-8, 
    alpha.0 = 1, 
    alpha.scheduler = NULL
 ){
  
  history <- data.frame(
    iter = numeric(0), 
    f = numeric(0),
    conv = numeric(0),
    alpha = numeric(0)
  )

  # If `alpha.scheduler` is NULL, use a constant value.
  if (is.null(alpha.scheduler)) {
    alpha.scheduler <- function(alpha, val, newval) alpha
  }
  
  # Define `step` function.
  step <- function(x, alpha) x + (alpha * g(x))
  
  # Perform a initial step and save it.
  i <- 1
  x <- x0  
  alpha <- alpha.0 
  newx <- step(x, alpha) 
  conv <- sum((x - newx) * (x - newx))
  history[nrow(history) + 1, ] <- c(i, f(newx), conv, alpha)
  
  while(conv > eps){
    
    # Update `alpha`.
    alpha <- alpha.scheduler(alpha, f(x), f(newx))
    
    # Update `x` to previous iteration's new value.
    x <- newx
    newx <- step(x, alpha)
    
    # Update convergence criterion.
    conv <- sum((x-newx)*(x-newx))
    
    # Save current iteration.
    i <- i + 1
    history[nrow(history) + 1, ] <- c(i, f(newx), conv, alpha)
    
  }
  
  return(list(x=newx, history=history))
  
}