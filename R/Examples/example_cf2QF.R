# Example (Generate specified quantiles from distribution specified by CF)
cf <- function(t) exp(-t^2/2)
prob <- matrix(runif(200), nrow = 100, ncol = 2)
result <- cf2QF(cf = cf, prob = prob)
qf <- result$qf
prob <- result$pr

# EXAMPLE (Generate Random Sample from distribution specified by CF)
  cf <- function(t) cf_Logistic(t)
  n  <- 10000
  p  <- runif(n)
  r  <- cf2QF(cf,p)
  hist(r,100)

# EXAMPLE (ECDF from distribution specified by CF)
  cf  <- function(t) cf_Logistic(t)
  n   <- 10000
  p   <- runif(n)
  r   <- cf2QF(cf,p)$qf
  pr  <- cf2QF(cf,p)$pr
  cdf <- cf2QF(cf,p)$cdf
  x   <- cf2QF(cf,p)$x
 plot(x,cdf)
