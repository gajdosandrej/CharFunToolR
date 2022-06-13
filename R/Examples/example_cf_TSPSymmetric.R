## EXAMPLE 1
# CF of the symmetric TSP distribution with theta =3/2 on (-1,1)
theta <- 3/2
t <- seq(from = -50,
         to = 50,
         length.out =501)
plotReIm(function(t)
        cf_TSPSymmetric(t, theta),
        t,
        title = "CF of the symmetric TSP distribution  on (-1,1)")



##EXAMPLE2
# PDF/CDF of the symmetric TSP distribution on (-1,1)
thet <- 3/2
cf <- function(t)
        cf_TSPSymmetric(t, theta)
x <- seq(-1,1,length.out = 101)
xRange <- 2
options <- list()
options.N <-2^8
options.dt <- 2*pi/xRange
result <- cf2DistGP(cf, x, c(), options)

##EXAMPLE 3
# CF of the weighted linear combination of TSP RVs
theta <- c(1, 2, 3, 4, 5)/2
mu <- c(1, 2, 0, 0, 0)
sigma <- c(1, 2, 3, 4, 5)/5
coef <- 1/5
t <- seq(-50,50, length.out=  501)
plotReIm(function(t)
        cf_TSPSymmetric(t, theta, mu, sigma, coef),
        t,
        title = "CF of the weighted linear combination of TSP RVs")

## EXAMPLE 4
# CDF/PDF of the weighted linear combination of TSP RVs
thet <- c(1, 2, 3, 4, 5)/2
mu <- 0
sigma <- c(5, 4, 3, 2, 1)
coef <- 1/5
t <-seq(-50, 50, length.out = 501)
cf <- function(t)
        cf_TSPSymmetric(t, theta, mu, sigma, coef)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options.N = 2^12
result <- cf2DistGP(cf, c(), prob, options)
