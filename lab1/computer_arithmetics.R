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