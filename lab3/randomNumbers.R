# TASK A
# ------

triangleDensity <- function(x) {
  ifelse(
    (x < -1) | (x > 1),
    0,
    ifelse((-1 <= x) & (x <= 0), x + 1, 1 - x)
  )
}


# The value of `c` by which we scale the uniform density should be as small
# as possible to minimize the frequency with which we reject the candidate 
# points. This requires determination of the maximum value of the triangle 
# density, which we can compute easily by evaluating at 0.
c <- triangleDensity(0)

# To generate the deviates from the triangle using the uniform majorizing 
# density (i.e., the envelope), we write the following statements.
y <- runif(10000, -1 , 1)  # Generate y from the majorizing (envelope) distribution.
u <- runif(10000)          # Generate u from a uniform (0, 1) distribution.
x <- na.omit(ifelse(u <= triangleDensity(y), y, NA))


# TASK B
# ------

# Let Y be a random variable following the triangle distribution from lecture 3
# with CDF F and F^-1 = 1 - sqrt(1 - y). Then, a sample from this distribution
# can be obtained as follows.
rtriang <- function(n, scale = 1) (1 - sqrt(1 - runif(n))) * scale

# Next, we can create a new random variable X using composition sampling
# based on Y and -Y.
n <- 10000
prob <- c(.3, .7)
scalers <- sample(c(1, -1), size = n, replace = TRUE, prob = prob)
x <- rtriang(n, scalers)


# TASK C
# ------

n <- 10000
u1 <- runif(n)
u2 <- runif(n)
x <- u1 - u2