# EXAMPLE 1
cf <- function(t) exp(-t^2/2)
options<-list()
options$isPlot <- FALSE
options$isInterp <- TRUE
result   <- cf2DistGP(cf,c(),c(),options)
xGiven   <- result$x
cdfGiven <- result$cdf
cdf <- function(x) InterpCDF(x,xGiven,cdfGiven)
x   <- seq(-10,10,length.out=1001)
plot (x,cdf(x),"l")

# EXAMPLE 2
cf <- function(t) exp(-t^2/2)
options <- list()
options$isPlot <- FALSE
options$isInterp <- TRUE
result <- cf2DistGP(cf = cf, options = options)
cdf    <- function(x) InterpCDF(x,result)
x      <- seq(-10,10,length.out=1001)
plot(x,cdf(x),"l")
