# EXAMPLE 1:
#CF of the log-transformed LRT statistic for testing compound symmetry
t  <- seq(-2,2,length.out=201)
n  <- 10
p  <- 5
type <- 'standard'
plotReIm(function(t) cfTest_CompoundSymmetry(t,n,p,type),t,
title='CF of test statistic for testing compound symmetry')


# EXAMPLE 2:
#CF of the log-transformed modified LRT statistic for compound symmetry
t  <- seq(-10,10,length.out=201)
n  <- 10
p  <- 5
type <- 'modified'
plotReIm(function (t) cfTest_CompoundSymmetry(t,n,p,type),t,
title='CF of test statistic for testing compound symmetry')


# EXAMPLE 3:
#PDF/CDFF/QF of the log-transformed LRT for testing compound symmetry
n  <- 10
p  <- 5
type <- 'standard'
cf   <- function(t) cfTest_CompoundSymmetry(t,n,p,type)
x    <- seq(0,25,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)

# EXAMPLE 4:
#PDF/CDFF/QF of the log-transformed modified LRT for compound symmetry
n  <- 10
p  <- 5
type <- 'modified'
cf   <- function(t) cfTest_CompoundSymmetry(t,n,p,type)
x    <- seq(0,5,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
