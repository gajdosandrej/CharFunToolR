## EXAMPLE 1
# CF of a linear combination of K=100 independent IGamma RVs
coef <- 1/(((1:50) - 0.5)*pi)^2
plot(1:50, coef, xlab="", ylab="", type = "p", pch = 20, col = "blue", cex = 1, main = expression('Coefficients of the linear combination of' ~ chi^2 ~ 'RVs with DF=1'))
lines(1:50, coef, col = "blue")
alpha <- 5/2
beta <- 2
t <- seq(from = -100,to = 100,length.out = 201)
plotGraf(function(t)
  cf_InverseGamma(t,alpha,beta,coef), t,
  title = "Characteristic function of the linear combination of IGamma RVs")

## EXAMPLE 2
# PDF/CDF from the CF by cf2DistGP
alpha <- 5/2
beta <- 1/2
coef <- 1/(((1:50) - 0.5)*pi)^2
cf <- function(t) cf_InverseGamma(t,alpha,beta,coef)
options <- list()
options$N <- 2^10
options$xMin = 0
x <- seq(from = 0,to = 4,length.out = 201)
prob <- c(0.9, 0.95, 0.99)
result <- cf2DistGP(cf,x,prob,options)
result
