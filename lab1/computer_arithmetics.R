myvar <- function(x){
  
  max_x <- max(x)
  x <- x / max_x
  
  
  n <- length(x)
  squaresum <- sum(x^2)
  sum_x <- sum(x)
  res <- (squaresum - (sum_x^2)/n)/(n-1)
  
  
  return(res * (max_x)^2)
}



myvar_res <- numeric(100)
var_res <- numeric(100)
Y <- numeric(100)
i_values <- 1:100
set.seed(12345)
for (i in i_values) {
  x <- rnorm(10000,  10^8, 1)
  myvar_res[i] <- myvar(x)
  var_res[i] <- var(x)
  Y[i] <- myvar_res[i] - var_res[i]
  
}

plot(i_values,Y,type = "l", main = "Difference between myvar(Xi) and var(Xi)",
     xlab = "i (Subset)", ylab = "Difference")


myvar_stable <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sum_sq_diff <- sum((x - mean_x)^2)
  return(sum_sq_diff / (n - 1))
}


myvar_stable_res <- numeric(100)
var_res_2 <- numeric(100)
Y_stable <- numeric(100)
for (i in i_values) {
  x <- rnorm(10000,  10^8, 1)
  myvar_stable_res[i] <- myvar_stable(x)
  var_res_2[i] <- var(x)
  Y_stable[i] <- myvar_stable_res[i] - var_res_2[i]
  
}

plot(i_values,Y_stable,type = "l", main = "Difference between myvar_stable(Xi) and var(Xi)",
     xlab = "i (Subset)", ylab = "Difference",ylim = c(-0.01,0.01))

# When the mean of the data is very large (e.g.,10^8) but the variance is small (e.g., 1), 
# the variance formula used in myvar function involves subtracting two very large numbers.
# This subtraction can cause catastrophic cancellation, where significant digits are lost due to floating-point precision limits.
