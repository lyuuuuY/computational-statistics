g_fun <- function(x){
  numerator <- log(x+1)
  denominator <- x^(3/2) + 1
  res <- numerator/denominator
  return(res)
}

x <- seq(0, 4, length.out = 100)
plot(x,g_fun(x),type = "b")

#part b

deriv_g <- function(x){
  u <- log(x+1)
  v <- x^(3/2) + 1
  deriv_u <- 1/(x+1)
  deriv_v <- (3/2) * x^(1/2)
  numerator <- deriv_u*v - u*deriv_v
  denominator <- v^2
  res <- numerator/denominator
  return(res)
}

plot(x,deriv_g(x),type = "b")
abline(h = 0, col = "red", lwd = 2, lty = 2)

#bisection

bisection <- function(g, g.prime, a, b, tol = 0.001) {
  xt <- c(0,0)
  xt[1] <- (a+b)/2
  while (abs(xt[2] - xt[1])>tol) {
    
    print(sprintf("x: %s", xt[1]))
    print(sprintf("g(x): %s", g(xt[1])))
    
    deriv_a <- g.prime(a)
    deriv_b <- g.prime(b)
    deriv_xt <- g.prime(xt[1])
    if((deriv_a*deriv_xt)<0){
      b <- xt[1]
      xt[2] <- xt[1]
      xt[1] <- (a+b)/2
    }
    else{
      a <- xt[1]
      xt[2] <- xt[1]
      xt[1] <- (a+b)/2
    }
  }
  x_value <- xt[1]
  max_value <- g(x_value)
  
}


# secant


secant <- function(g, g.prime, initial.guess, tol = 0.001) {
  xt_secant <- c(initial.guess,0)
  while (abs(xt_secant[2] - xt_secant[1]) > tol){
    
    print(sprintf("x: %s", xt_secant[1]))
    print(sprintf("g(x): %s", g(xt_secant[1])))
    
    deriv_xt <- g.prime(xt_secant[1])
    deriv_xt_minus1 <- g.prime(xt_secant[2])
    x_tplus1 <- xt_secant[1] - deriv_xt * (xt_secant[1]-xt_secant[2])/(deriv_xt-deriv_xt_minus1)
    xt_secant[2] <- xt_secant[1]
    xt_secant[1] <- x_tplus1
  }
  x_value_secant <- xt_secant[1]
  max_value_secant <- g(x_value_secant)  
}












