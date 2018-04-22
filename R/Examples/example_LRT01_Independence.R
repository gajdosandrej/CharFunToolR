## EXAMPLE
# LRT for testing hypothesis about independence
# Null distribution of the minus log-transformed LRT statistic
n <- 30                  # sample size
p <- c(3,4,5,6,7)        # dimensions of X_k,k = 1,...,m where m = 5
# W <- vector()          # observed value of W = -log(Lambda)
options <- list()
# options$coef <- -1
options$prob <- c(0.9,0.95,0.99)
output <- LRT01_Independence(n = n, p = p, options = options)
str(output)
