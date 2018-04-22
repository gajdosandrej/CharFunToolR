## EXAMPLE
# LRT for testing hypothesis on sphericity of covariance matrix
# Null distribution of the minus log-transformed LRT statistic
n <- 30                    # total sample size
p <- 8                     # dimension of X ~ N_p(mu,Sigma)
# W <- vector()            # observed value of W = -log(Lambda)
options <- list()
# options$coef = -1
options$prob <- c(0.9, 0.95, 0.99)
output <- LRT05_Sphericity(n = n, p = p, options = options)
str(output)
