## EXAMPLE 1
# PMF/CDF of the binomial RVV specified bz its CF
 n<-10
 p<-0.25
 cf<-function(t) cfN_Binomial(t,n,p)
  xMin<- 0
  xMax<-n
  xDelta<-1
  options$isPlot<-TRUE
  result<-cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
  
  ## EXAMPLE 2   
  # PMF/CDF of the convolved discrete RV specified by its CF
  n<-10
  p<-0.25
  cf_Bino<-function(t) cfN_Binomial(t,n,p)
  N<-5
  cf<-function(t) cf_Bino(t)^N
  xMin<-0
  xMax<-N*n
  xDelta<-1
  options$isPlot=TRUE
  result<-cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
  
  ## EXAMPLE 3
  # PMF/CDF of the mean of IID RV with Poisson distribution
  lambda<-5
  cf_Pois<-function(t) cfN_Poisson(t,lambda)
  N<-10
  cf<-function(t) cf_Pois(t/N)^N
  xMin<-0
  xMax<-10
  xDelta<-1/N
  options$isPlot=TRUE
  result<-cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
  
  ## EXAMPLE 4
  # PMF/CDF of the convolved discrete RV specified by its CF
  # Here we consider convolutions of discrete RV defined on{0,1,2}
  # with probabilities p=[0.2,0.5,0.3]
  supp<-c(0,1,2)
  prob<-c(0.2,0.5,0.3)
  cf_X<-function(t) cfE_DiracMixture(t,supp,prob)
  N<-5
  cf<-function(t) cf_X(t)^N
  xMin<-0
  xMax<-max(supp)*N
  xDelta<-1
  options$isPlot<-TRUE
  result<-cf2PMF_FFT(cf,xMin,xMax,xDelta,options)
  
  ## EXAMPLE 5
  # PMF/CDF of the exact bootrap mean distribution specified by its CF
  data<-c(-2,0,0,0,0,0,0,0,1,4)
  N<-length(data)
  cf<- function(t) cfE_Empirical(t/N,data)^N
  xMin<-min(data)
  xMax<-max(data)
  xDelta<-1/N
  options$isPlot<-TRUE
  result<-cf2PMF_FFT(cf,xMin, xMax, xDelta, options)
  
  