# EXAMPLE 1:
#CF of the log-transformed LRT statistic for equal covariances
t  <- seq(-1,1,length.out=201)
n  <- 10
p  <- 5
q  <- 3
type <- 'standard'
plotReIm(function (t) cfTest_EqualityCovariances(t,n,p,q,type),t,
title='CF of test statistic for testing equal covariances',ylab='cf(t)')

# EXAMPLE 2:
#CF of the log-transformed modified LRT statistic for equal covariances
t  <- seq(-5,5,length.out=201)
n  <- 10
p  <- 5
q  <- 3
type <- 'modified'
plotReIm(function(t) cfTest_EqualityCovariances(t,n,p,q,type),t,
title='CF of test statistic for testing equal covariances')

#EXAMPLE 3:
#PDF/CDFF/QF of the log-transformed LRT for equal covariances
n    <- 10
p    <- 5
q    <- 3
type <- 'standard'
cf   <- function(t) cfTest_EqualityCovariances(t,n,p,q,type)
x    <- seq(0,50,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options<-list()
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)

# EXAMPLE 4:
#PDF/CDFF/QF of the log-transformed modified LRT for equal covariances
n    <- 10
p    <- 5
q    <- 3
type <- 'modified'
cf   <- function(t) cfTest_EqualityCovariances(t,n,p,q,type)
x    <- seq(0,10,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
