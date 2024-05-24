##   EXAMPLE 1:
##Empirical CF of the sample mean
    data <- c(1, 1, 2, 2, 2, 3, 3, 9, 105, 105, 106, 106, 106, 107, 107)
    result <- cfE_SampleMean(c(),data)
    x  <- unique(result$X)
    cf  <- result$cf
    t   <- seq(-0.5,0.5,length.out=1001)
    plotReIm(function(t) cf(t),ylab = "cf(t)",
             t, title = 'Empirical CF of the Sample Mean')

##   EXAMPLE 2:
## Empirical CF and PMF/CDF of the sample mean
    data   <- c(-2, -1, 0, 1, 2, 3, 4)
    counts <- c( 2,  2, 3, 1, 2, 3, 2)
    resultCF <- cfE_SampleMean(c(),data,counts)
    x   <- unique(resultCF$X)
    cf  <- resultCF$cf
    t   <- seq(-100,100,length.out=1001)

    plotReIm(function(t) cf(t),
             t, title = 'Empirical CF of the Sample Mean')
    n <- resultCF$n
    delta <-1/n
    resultFFT <- cf2PMF_FFT(cf,-2,4,delta)
 ## Numerical inversion of the smoothed CF by cf2DistGP
    cfSmooth <- function(t) cf(t) * cf_Normal(t*delta)
    resultGP  <- cf2DistGP(cfSmooth)

##   EXAMPLE 3:
   ## Empirical distribution of the difference of the sample means
## Test of the hypothesis that the population means are equal
    dataA   <- c(1, 1, 2, 2, 2, 3, 3, 9, 105, 105, 106, 106, 106, 107, 107)
    dataB   <- c(5, 5, 6, 6, 6, 7, 7, 99, 101, 101, 102, 102, 102, 103, 103)
    cfA     <- function(t) cfE_SampleMean(t,dataA)
    cfB     <- function(t) cfE_SampleMean(t,dataB)
    cfDiff  <- function(t) cfA(t)$cf* cfB(-t)$cf
    minDiff <- min(dataA) - max(dataB)
    maxDiff <- max(dataA) - min(dataB)
    delta   <- 1/15
    result <- cf2PMF_FFT(cfDiff,minDiff,maxDiff,delta)

##   EXAMPLE 4:
## Empirical distribution of the difference of the sample means
# Test of the hypothesis that the population means are equal
    dataA   <- c(1, 1, 2, 2, 2, 3, 3, 9, 105, 105, 106, 106, 106, 107, 107)
    dataB  <- c(5, 5, 6, 6, 6, 7, 7, 99, 101, 101, 102, 102, 102, 103, 103)
    options<-list()
    options$isOgive <- TRUE
    cfA     <- function(t) cfE_SampleMean(t,dataA,c(),options)
    cfB     <- function(t) cfE_SampleMean(t,dataB,c(),options)
    cfDiff  = function(t) cfA(t)$cf * cfB(-t)$cf
    result <- cf2DistGP(cfDiff,c(),c(0.25,0.975))
