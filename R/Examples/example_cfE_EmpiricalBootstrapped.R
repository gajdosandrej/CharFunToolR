#EXAMPLE 1 (Bootstrapped empirical CF based on given data)
set.seed(101)
n <- 1000
Data <- c(rnorm(3 * n, 5, 0.2), rt(n, 3),rchisq(n, df=1))

  dataN <- length(Data)
  randID <-pracma:: randi(dataN,dataN,1);
  t <- seq(-50, 50, length.out = 2^10)
  plotReIm(function(t) cfE_EmpiricalBootstrapped(t, Data,t(randID)),
           t, title = 'Bootstrapped empirical CF')


# EXAMPLE 2 (PDF/CDF of the compound bootstrapped empirical distributions)
  set.seed(101)
  lambda <- 25
  nN <- 100
  Ndata <- rpois(nN,lambda)
  mu <- 0.1
  sigma <- 2
  nX <- 15000
  Xdata <- rlnorm(nX,mu,sigma)
  randID <- vector()
  t   <-seq(-0.2,0.2,length.out = 2^10)
  cfX <-function(t) cfE_EmpiricalBootstrapped(t,Xdata,randID)
cf<-function(t) cfE_EmpiricalBootstrapped(t,Ndata,randID,cfX)

  plotReIm(function(t) cfE_EmpiricalBootstrapped(t,Ndata,randID,cfX),
           t, title = 'Bootstrapped empirical CF')
  x <- seq(0,1000,length.out = 501)
  prob <- c(0.9, 0.95)
  options <- list()
  options$N <- 2^12
  options$xMin <- 0
  options$SixSigmaRule <- 10
  result <- cf2DistGP(cf,x,prob,options)

# EXAMPLE 3 (PDF/CDF of the Stress-Strength reliability R = Pr(X<Y))
  set.seed(101)
  mu <- 0
  sigma <- 1
  n <- 50
  X <- rlnorm(n,mu,sigma)
  mu <- 2
  sigma <- 2
  n <- 20
  Y <- rlnorm(n,mu,sigma)
  nBoot <- 1000
  xCrit <- 0
  R     <- rep(0,nBoot)

  for (i in 1:nBoot){
      cfX <-function(t) cfE_EmpiricalBootstrapped(t,X)
      cfY <- function(t) cfE_EmpiricalBootstrapped(t,Y)
      cf  <- function(t) cfX(t) * cfY(-t)
      M <- cf2DistGP(cf,xCrit)
      R[i]<-M$cdf
  }
  bandwidth <- 0.02
  cf_KERNEL <- function(t) exp(-(bandwidth*t)^2/2)
  cfR <- function(t) cfE_Empirical(t,R) * cf_KERNEL(t)
  x <- seq(0,1,length.out=100)
  prob <-c(0.025, 0.5, 0.95, 0.975)
  options<-list()
  options$xMin <- 0
  options$xMax <- 1
  options$SixSigmaRule <- 10
  result <- cf2DistGP(cfR,x,prob,options)
