## EXAMPLE 1
# CF of the log-transformed non-central Beta RV with delta = 1, coef = -1
alpha <- 1
beta  <- 3
delta <- 1
coef  <- -1
t <- seq(from = -10,
         to = 10,
         length.out = 201)
functions_to_plot <- list(function(t) cf_LogRV_BetaNC(t, alpha, beta, delta, coef, type = 1),
                          function(t) cf_LogRV_BetaNC(t, alpha, beta, delta, coef, type = 2))
plotReIm2(functions_to_plot, list(t, t),
          title = 'CF of minus log-transformed Type I and II Beta RV')


type <- 1
plotReIm(function(t)
        cf_LogRV_BetaNC(t, alpha, beta, delta, coef, type),
        t,
        title = 'CF of minus log-transformed Type I and II Beta RV')
par(new=TRUE)
type <- 2
plotReIm(function(t)
        cf_LogRV_BetaNC(t, alpha, beta, delta, coef, type),
        t,
        title = 'CF of minus log-transformed Type I and II Beta RV')

## EXAMPLE 2
# CDF/PDF of the minus log-transformed non-central Beta RV with delta = 1
alpha <- 1
beta  <- 3
delta <- 1
coef  <- -1
cf <- function(t)
        cf_LogRV_BetaNC(t, alpha, beta, delta, coef)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)

## EXAMPLE 3
# CDF/PDF of the linear combination of minus log-transformed non-central Beta RVs
alpha <- c(1, 2, 3)
beta  <- c(3, 4, 5)
delta <- c(0, 1, 2)
coef  <- c(-1,-1,-1) / 3
cf <- function(t)
        cf_LogRV_BetaNC(t, alpha, beta, delta, coef)
options <- list()
options$xMin <- 0
result <- cf2DistGP(cf = cf, options = options)
