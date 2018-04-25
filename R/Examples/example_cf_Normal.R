## EXAMPLE 1
# CF of a linear combination of K=100 independent Norma RVs
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
plot(
        1:50,
        coef,
        xlab = "",
        ylab = "",
        type = "p",
        pch = 20,
        col = "blue",
        cex = 1,
        main = expression('Coefficients of the linear combination of' ~ chi ^ 2 ~ 'RVs with DF=1')
)
lines(1:50, coef, col = "blue")
mu <- seq(from = -3,
          to = 3,
          length.out = 50)
sigma <- seq(from = 0.1,
             to = 1.5,
             length.out = 50)
t <- seq(from = -100,
         to = 100,
         length.out = 2001)
plotGraf(function(t)
        cf_Normal(t, mu, sigma, coef),
        t,
        title = "Characteristic function of the linear combination of Normal RVs")

## EXAMPLE 2
# PDF/CDF from the CF by cf2DistGP
mu <- seq(from = -3,
          to = 3,
          length.out = 50)
sigma <- seq(from = 0.1,
             to = 1.5,
             length.out = 50)
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
cf <- function(t)
        cf_Normal(t, mu, sigma, coef)
options <- list()
options$N <- 2 ^ 10
options$SixSigmaRule <- 8
prob <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.950, 0.99)
result <- cf2DistGP(cf = cf, prob = prob, options = options)
result
