## EXAMPLE1
# CF of the Student t-distribution with df = 3
df <- 3
t <- seq(-5, 5, length.out = 501)
plotGraf(function(t)
        cfS_Student(t, df), t, title = "CF of the Student t-distribution with df = 3")

## EXAMPLE2
# PDF/CDF of the Student t-distribution with df = 3
df <- 3
cf <- function(t)
        cfS_Student(t, df)
x <- seq(-8, 8, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2 ^ 12
options$SixSigmaRule <- 30
result <- cf2DistGP(cf, x, prob, options)
