## EXAMPLE
# LRT for testing hypothesis on equality of populations
# Null distribution of the minus log-transformed LRT statistic
n <- 30                  # total sample size
p <- 8                   # dimension of X_k, k = 1,...,q where q = 5
q <- 5                   # number of populations
# W <- vector()          # observed value of W = -log(Lambda)
options <- list()
# options.coef = -1;
options$prob <- c(0.9, 0.95, 0.99)
output <- LRT04_EqualityPopulations(n = n, p = p, q = q, options = options)
str(output)
