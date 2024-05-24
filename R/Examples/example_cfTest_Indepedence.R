# EXAMPLE 1:
#CF of the log-transformed LRT statistic for indepedence
t  <- seq(-1,1,length.out=201)
n  <- 20
p  <- c(4,5,6)
q  <- c()
type <- 'standard'
plotReIm (function(t) cfTest_Indepedence(t,n,p,q,type),t,
          title='CF of LRT for testing indepedence of q normal populations')

# EXAMPLE 2:
#CF of the log-transformed LRT statistic  for indepedence
t  <- seq(-1,1,length.out=201)
n  <- 20
p  <- 5
q  <- 3
type <- 'standard'
plotReIm(function(t) ccfTest_Indepedence(t,n,p,q,type),t,
         title='CF of LRT for testing indepedence of q normal populations')

# EXAMPLE 3:
#CF of the log-transformed modified LRT statistic for independence
t  <- seq(-5,5,length.out=201)
n  <- 20
p  <- 5
q  <- 3
type <- 'modified'
plotReIm (function(t) cfTest_Indepedence(t,n,p,q,type),t,
title='CF of LRT for testing independence of q normal populations')

#EXAMPLE 4:
#PDF/CDFF/QF of the log-transformed LRT for independence
n    <-  20
p    <- 5
q    <- 3
type <- 'standard'
cf   <- function(t) cfTest_Indepedence(t,n,p,q,type)
x    <- seq(30,120,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)

#EXAMPLE 5:
#PDF/CDFF/QF of the log-transformed modified LRT for indepedence
n  <- 20
p    <- 5
q    <- 3
cf   <- function(t) cfTest_Indepedence(t,n,p,q)
x    <- seq(2,12,length.out=201)
prob <- c(0.9, 0.95, 0.99)
options$xMin <- 0
result <- cf2DistGP(cf,x,prob,options)
