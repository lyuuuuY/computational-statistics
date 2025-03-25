# --- ASSIGNMENT 1 ---
# We consider a case when a Gumbel distributed observation x is made. The 
# Gumbel distribution has a known scale parameter 1; the location parameter μ 
# is unknown.

# Functions definition.
lik <- function(mu, x) exp(-x + mu  - exp(-x + mu))
lik_derivative <- function(mu, x) lik(mu,x) * (1 - exp(-x + mu))

prior <- function(mu) exp(-mu^2 / 2)
prior_derivative <- function(mu) -mu * prior(mu)

posterior <- function(mu, x) lik(mu, x) * prior(mu)
posterior_derivative <- function(mu, x) lik(mu, x) * prior_derivative(mu) + lik_derivative(mu, x) * prior(mu)


# Observed value.
x <- 1.65

# >>> Question 1.1
# Plot the posterior density.
mu_values <- seq(-5, 5, length.out = 50)
posterior_values <- posterior(mu_values, x = x)
plot(mu_values, posterior_values)

# Plot also the derivative.
posterior_derivative_values <- posterior_derivative(mu_values, x = x)
plot(mu_values, posterior_derivative_values, type = "l")


# >>> Question 1.2
# Maximize the posterior using the secant method
mysecant <- function(a, b, FUN, tol = 1e-6) {
  
  step <- function(x_t_minus_1, x_t_minus_2) x_t_minus_1 - (FUN(x_t_minus_1) * (x_t_minus_1 - x_t_minus_2) / (FUN(x_t_minus_1) - FUN(x_t_minus_2)))

  x_t_minus_2 <- a
  x_t_minus_1 <- b
  x_t <- step(x_t_minus_1, x_t_minus_2)
  
  conv <- abs(x_t - x_t_minus_1)
  i <- 0
  
  while (conv > tol) {
    x_t_minus_2 <- x_t_minus_1
    x_t_minus_1 <- x_t
    x_t <- step(x_t_minus_1, x_t_minus_2)
    conv <- abs(x_t - x_t_minus_1)
    i <- i + 1
  }
  
  return(x_t)
}

objective <- function(x) posterior_derivative(x, 1.65)

# Choose two pairs of starting values: One pair which leads to convergence
# to the right global maximum at ˆμ and another pair not leading to convergence to that maximum.
# Explain what happens in the unsuccessful case.
success <- mysecant(-1, 1, objective)
failure <- mysecant(-2, -1, objective)


# >>> Question 1.3
make_metropolis_hasting_sampler <- function(x0, rproposed, dproposed, dtarget) {
  
  metropolis_hasting_sampler <- function(n) {
    
    x <- x0
    sample <- numeric(n)
    count <- 0
    
    for (i in 1:n) {
      
      # Generate y from the Y_t+1 | Y_t = x
      y <- rproposed(x)
      
      # Compute ratio
      r <- dtarget(y) * (dproposed(x, y) / (dtarget(x) * dproposed(y, x)))
      
      # Select or Reject sample y
      x <- ifelse(
        r >= 1, 
        {count <- count + 1; y},
        ifelse(runif(1) < r,  {count <- count + 1; y}, x))
      
      sample[i] <- x
    }
    return(sample)
  }
  
  return(metropolis_hasting_sampler)
}

rproposed <- function(x, a) x + runif(1, -a, a)
dproposed <- function(x, y, a) dunif(x, y - a, y + a)
dtarget <- function(x) posterior(x, 1.65)
n <- 1000


# a = 0.15
sampler_a015 <- make_metropolis_hasting_sampler(
  x0 = 0,
  rproposed = function(x) rproposed(x, a = 0.15),
  dproposed = function(x, y) dproposed(x, y, a = 0.15),
  dtarget = dtarget
)
sample_values_a015 <- sampler_a015(n)

# a = 1.5
sampler_a15 <- make_metropolis_hasting_sampler(
  x0 = 0,
  rproposed = function(x) rproposed(x, a = 1.5),
  dproposed = function(x, y) dproposed(x, y, a = 1.5),
  dtarget = dtarget
)
sample_values_a15 <- sampler_a15(n)

plot(sample_values_a15, type = "l")


# >>> Question 1.4
pr_mu_greater_than_zero <- sum(sample_values_a15 >= 0) / n




# --- ASSIGNMENT 2 ---

# >>> Question 2.1
# Visualize the function to maximize and mark the global maximum.
f <- function(x) (x^2 / exp(x)) - 2 * exp(-(9 * sin(x)) / (x^2 + x + 1))
res <- optimize(f, maximum = TRUE, interval = c(0, 5))
x <- seq(0, 10, length.out = 100)
plot(x, f(x), type = "l")
points(res$maximum, f(res$maximum), col="red")


# >>> Question 2.2









