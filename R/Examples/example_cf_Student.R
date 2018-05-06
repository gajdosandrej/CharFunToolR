## EXAMPLE 1
# CF of a linear combination of independent Student's t RVs
coef <- 1 / (1:50)
df   <- seq(50, 1, -1)
t <- seq(from = -1,
         to = 1,
         length.out = 201)
plotReIm(function(t)
        cf_Student(t, df, coef),
        t,
        title = "Characteristic function of the linear combination of t RVs")

## EXAMPLE 2
# CDF/PDF of a linear combination of independent Student's t RVs
coef <- 1 / (1:50)
df <- seq(50, 1, -1)
cf <- function(t)
        cf_Student(t, df, coef)
x <- seq(from = -50,
         to = 50,
         length.out = 100)
prob <- c(0.9, 0.95, 0.975, 0.99)
options <- list()
options$N <- 2 ^ 12
result <- cf2DistGP(cf, x, prob, options)
result
