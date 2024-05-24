
##EXAMPLE 1:
## Empirical CF and PMF/CDF of the sample median
    dataA   <- c(1, 1, 2, 2, 2, 3, 3, 9, 105, 105, 106, 106, 106, 107, 107)
    result <- cfE_SampleMedian(c(),dataA)
    x   <- result$Xunique
    pmf <- result$PMFmedian
    cdf <- result$CDFmedian
    cf  <- result$cf
    t   <- seq(-10,10,length.out=1001)
    plot(x,pmf,"h",main='Empirical PMF of the Sample Median')
    plot(x,cdf,"s",main='Empirical CDF of the Sample Median')
       plotReIm(function(t) cf(t),title='Empirical CF of the Sample Median')
##   EXAMPLE 2:
##Empirical distribution of the difference of the sample medians
# # Test of the hypothesis that the population medians are equal
    dataA   <-c(1, 1, 2, 2, 2, 3, 3, 9, 105, 105, 106, 106, 106, 107, 107)
    dataB   <-c(5, 5, 6, 6, 6, 7, 7, 99, 101, 101, 102, 102, 102, 103, 103)
    cfA     <- function(t) cfE_SampleMedian(t,dataA)
    cfB     <- function(t) cfE_SampleMedian(t,dataB)
    cfDiff  <- function(t) cfA(t)$cf * cfB(-t)$cf
    minDiff <- min(dataA) - max(dataB)
    maxDiff <- max(dataA) - min(dataB)
    delta   <- 1
    result <- cf2PMF_FFT(cfDiff,minDiff,maxDiff,delta)
