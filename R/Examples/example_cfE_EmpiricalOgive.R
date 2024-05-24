## EXAMPLE1 (Ogive ECF - a weighted mixture of independent Uniform RVs)
set.seed(101)
n <- 1000
Data <- c(rnorm(3 * n, 5, 0.2), rt(n, 10))
Edges <- -6:6
Counts <- hist(x = Data, breaks = Edges)$counts
Xbins <- Edges
Xfreq <- Counts / sum(Counts)
t <- seq(-10, 10, length.out = 501)
bins <- Xbins
freq <- Xfreq
plotReIm(function(t) cfE_EmpiricalOgive(t, Xbins, Xfreq),
         t, title = 'Ogive ECF of the grouped (histogram) data')

## EXAMPLE2 (PDF/CDF of the Ogive distribution)
set.seed(101)
n <- 1000
Data <- c(rnorm(3 * n, 5, 0.2), rt(n, 10))
Edges <- -6:6
Counts <- hist(x = Data, breaks = Edges)$counts
Xbins <- Edges
Xfreq <- Counts / sum(Counts)
cf <- function(t) cfE_EmpiricalOgive(t, Xbins, Xfreq)
x <- seq(-5, 7, length.out = 101)
prob <- c(0.9, 0.95)
options <- list()
options$N <- 2^12
options$SixSigmaRule <- 8
result <- cf2DistGP(cf, x, prob, options)

## EXAMPLE3 (PDF/CDF of the compound Poisson-Ogive distribution)
set.seed(101)
n <- 1000
Data <- c(rnorm(3 * n, 5, 0.2), rt(n, 10))
Edges <- -6:6
Counts <- hist(x = Data, breaks = Edges)$counts
Xbins <- Edges
Xfreq <- Counts / sum(Counts)
cfX <- function(t) cfE_EmpiricalOgive(t, Xbins, Xfreq)
lambda <- 15
cf <- function(t) cfN_Poisson(t, lambda, cfX)
x <- seq(0, 100, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x, prob, options)
