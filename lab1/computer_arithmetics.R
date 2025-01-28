myvar <- function(x){
  
  max_x <- max(x)
  x <- x / max_x
  
  
  n <- length(x)
  squaresum <- sum(x^2)
  sum_x <- sum(x)
  res <- (squaresum - (sum_x^2)/n)/(n-1)
  
  
  return(res * (max_x)^2)
}



myvar_res <- numeric(10000)
var_res <- numeric(10000)
Y <- numeric(10000)
i_values <- 1:10000
set.seed(12345)
x <- rnorm(10000,  10^8, 1)
for (i in i_values) {
  sub_x <- x[1:i]
  myvar_res[i] <- myvar(sub_x)
  var_res[i] <- var(sub_x)
  Y[i] <- myvar_res[i] - var_res[i]
  
}

plot(i_values,Y, main = "Difference between myvar(Xi) and var(Xi)",
     xlab = "i (Subset)", ylab = "Difference")



plot(1:100, Y[1:100], type = "b", col = "blue",
     main = "Details of the differences in the first 100 samples",
     xlab = "i", ylab = "Y_i")

myvar_stable <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sum_sq_diff <- sum((x - mean_x)^2)
  return(sum_sq_diff / (n - 1))
}


myvar_stable_res <- numeric(10000)
var_res_2 <- numeric(10000)
Y_stable <- numeric(10000)
for (i in i_values) {
  sub_x <- x[1:i]
  myvar_stable_res[i] <- myvar_stable(sub_x)
  var_res_2[i] <- var(sub_x)
  Y_stable[i] <- myvar_stable_res[i] - var_res_2[i]
  
}

plot(i_values,Y_stable, main = "Difference between myvar_stable(Xi) and var(Xi)",
     xlab = "i (Subset)", ylab = "Difference",ylim = c(-0.01,0.01))

# When the mean of the data is very large (e.g.,10^8) but the variance is small (e.g., 1), 
# the variance formula used in myvar function involves subtracting two very large numbers.
# This subtraction can cause catastrophic cancellation, where significant digits are lost due to floating-point precision limits.
