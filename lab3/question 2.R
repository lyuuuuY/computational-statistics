mu <- c(0,0)
sigma <- matrix(c(0.6,0,0,0.6),nrow = 2)

box_muller <- function(n,mu,sigma){
  u <- runif(n,0,1)
  A <- runif(n,0,2*pi)
  x1 <- sqrt(-2*log(u)) * cos(A)
  x2 <- sqrt(-2*log(u)) * sin(A)
  x_standard <- cbind(x1,x2)
  w <- chol(sigma)
  Y <- x_standard%*%t(w) + matrix(rep(mu, each = n), ncol = 2)
  return(Y)
}

system.time(box_muller(10000000,mu,sigma))


library(mvtnorm)
system.time(rmvnorm(10000000,mu,sigma,method = "chol"))


mix_normal <- function(n,prob,mu1,mu2,sigma1,sigma2){
  
  # get the number of observations of each density
  density1_n <- floor(n*prob)
  density2_n <- n-density1_n
  
  sample1 <- box_muller(density1_n,mu1,sigma1)
  sample2 <- box_muller(density2_n,mu2,sigma2)
  
  # add lable to each component
  sample1 <- cbind(sample1,1)
  sample2 <- cbind(sample2,2)
  
  samples <- rbind(sample1,sample2)
  return(samples)
}

mu1 <- c(0,0)
mu2 <- c(1.5,1.2)
sigma1 <- matrix(c(0.6,0,0,0.6),nrow = 2)
sigma2 <- matrix(c(0.5,0,0,0.5),nrow = 2)

new_samples <- mix_normal(1000,0.5,mu1,mu2,sigma1,sigma2)
sample1 <- new_samples[which(new_samples[,3] == 1), ]
sample2 <- new_samples[which(new_samples[,3] == 2), ]
plot(sample1[,1],sample1[,2],col="red",ylim = c(-3,3),xlim = c(-3,3))
points(sample2[,1],sample2[,2],col="blue")
