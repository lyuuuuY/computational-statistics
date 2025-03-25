em <- function(data, eps = 0.0001){
  n      <- length(data)
  pi     <- matrix(NA, ncol = 3, nrow = n)  
  
  # Random initialization of the parameters
  p <- rep(1/3, 3)
  sigma <- rep(sd(data)*2/3, 3)
  mu <- c(mean(data) - sigma[1]/2, mean(data), mean(data) + sigma[1]/2)
  pv <- c(p, mu, sigma)
  cc  <- eps + 100
  mu_list <- c(mu)
  sigma_list <- c(sigma)
  
  while (cc > eps){
    # save previous parameter vector
    pv1  <- pv                 
    
    # E step 
    for (j in 1:n){
      pi1 <- p[1]*dnorm(data[j], mean=mu[1], sd=sigma[1])
      pi2 <- p[2]*dnorm(data[j], mean=mu[2], sd=sigma[2])
      pi3 <- p[3]*dnorm(data[j], mean=mu[3], sd=sigma[3])
      
      total <- pi1+pi2+pi3
      pi[j,] <- c(pi1/total, pi2/total, pi3/total)
    }
    
    # M step #
    p <- colSums(pi)/n
    mu[1] <- sum(pi[, 1]*data)/(p[1]*n)
    mu[2] <- sum(pi[, 2]*data)/(p[2]*n)
    mu[3] <- sum(pi[, 3]*data)/(p[3]*n)
    sigma[1] <- sqrt(sum(pi[, 1]*((data-mu[1])^2))/(p[1]*n))
    sigma[2] <- sqrt(sum(pi[, 2]*((data-mu[2])^2))/(p[2]*n))
    sigma[3] <- sqrt(sum(pi[, 3]*((data-mu[3])^2))/(p[3]*n))
    
    # stopping criteria
    pv <- c(p, mu, sigma)
    cc     <- sum((pv - pv1)^2) / sum(pv1^2)  # a convergence criterion, maybe not the best one
    mu_list <- rbind(mu_list, mu)
    sigma_list <- rbind(sigma_list, sigma)
    
  }
  res <- list(final_probs = pv[1:3], 
              final_mus = pv[4:6],
              final_sigmas = pv[7:9],
              mu = mu_list, 
              sigma = sigma_list)
  return(res)
}


data <- c(0.1, 0.5, 0.7, 1.1, 2.5, 3.4, 3.5, 3.9, 4.0)
res2 <- em(data)
