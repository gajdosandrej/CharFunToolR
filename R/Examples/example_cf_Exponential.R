## EXAMPLE 1
# CF of the Exponential distribution with lambda = 5
lambda <- 5
t <- seq(from = -50,
         to = 50,
         length.out = 501)
plotGraf(
        function(t)
                cfX_Exponential(t, lambda),
        t,
        title = expression('CF of the Exponential distribution with' ~ lambda ~ '= 5')
)

## EXAMPLE 2
# CF of a linear combination of independent Exponential RVs
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
lambda <- 5
t <- seq(from = -100,
         to = 100,
         length.out = 201)
plotGraf(function(t)
        cfX_Exponential(t, lambda, coef),
        t,
        title = "CF of a linear combination of EXPONENTIAL RVs")

## EXAMPLE 3
# PDF/CDF of the compound Binomial-Exponential distribution
n <- 25
p <- 0.3
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
lambda <- 5
cfX <- function(t)
        cf_Exponential(t, lambda, coef)
cf <- function(t)
        cfN_Binomial(t, n, p, cfX)
x <- seq(from = 0,
         to = 5,
         length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
options$N <- 2 ^ 12
options$SixSigmaRule <- 15
result <- cf2DistGP(cf, x, prob, options)
