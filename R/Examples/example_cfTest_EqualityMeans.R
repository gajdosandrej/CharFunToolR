# EXAMPLE 1:
#CF of the log-transformed LRT statistic for equal means
t  <- seq(-2,2,length.out=201)
n  <- 12 # say n = 3+4+5
p  <- 5
q  <- 3
type <- 'standard'
plotReIm (function(t) cfTest_EqualityMeans(t,n,p,q,type),t,
title='CF of test statistic for testing equal means')

# EXAMPLE 2:
#CF of the log-transformed modified LRT statistic for equal means
t  <- seq(-10,10,length.out=201)
n  <- 12 # say n = 3+4+5
p  <- 5
q  <- 3
type <- 'modified'
plotReIm(function(t) cf = cfTest_EqualityMeans(t,n,p,q,type),t,
title='CF of test statistic for testing equal means')

#EXAMPLE 3:
#PDF/CDFF/QF of the log-transformed LRT for equal means
n  <- 12 # say n = 3+4+5
p    <- 5
q    <- 3
type <- 'standard'
cf   <- function(t) cfTest_EqualityMeans(t,n,p,q,type)
x    <- seq(0,30,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)

#EXAMPLE 4:
#PDF/CDFF/QF of the log-transformed modified LRT for equal means
n  <- 12 # say n = 3+4+5
p    <- 5
q    <- 3
type <- 'modified'
cf   <- function(t) cfTest_EqualityMeans(t,n,p,q,type)
x    <- seq(0,5,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
