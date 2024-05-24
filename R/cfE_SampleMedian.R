#' @title Empirical characteristic function of the sample median based on the
#'    observed data - the realization of the random sample
#'
#' @description
#'   The distribution of the sample median is estimated by its empirical
#'   distribution, which is defined as the distribution of the kth-order
#'   statistics derived from the observed empirical cumulative distribution
#'   function (ECDF). The order of the median is given by \eqn{k = (n+1)/2}.
#'
#'   Note that for the continuous distribution F, the kth-order statistic
#'   has a beta distribution and its CDF values at \eqn{x} are given as \code{CDF(x) = pbeta(F(x),k,n+1-k)}.
#'   Surprisingly, this is also true for discrete distributions, such as the ECDF,
#'   when evaluated at the support points,i.e. \code{CDF (x) = pbeta(ECDF(x),k,n+1-k)}.
#'   The empirical CDF of the sample median is then \code{CDF (x) = pbeta(ECDF(x),(n+1)/2,(n+1)/2)}.
#'
#'   When the sample size \eqn{n} is an odd integer, CDF corresponds to the
#'   bootstrap distribution of the sample median (i.e. the limiting
#'   bootstrap distribution calculated from the large number of bootstrap
#'   samples \eqn{B -> infinity}).
#'
#'   When the sample size \eqn{n} is an even integer, the sample median is usually
#'   defined as the mean of the \eqn{kth} and the \eqn{(k+1)th} order statistics, where
#'  \eqn{ k = n/2}. However, \code{CDF (x) = pbeta(ECDF(x),(n+1)/2,(n+1)/2)} is also
#'   well defined for \eqn{n} an even integer (but this CDF would be slightly
#'   different from the bootstrap distribution as it could be concentrated
#    also on points that are out of the observed data). Therefore, we define
#'   here the empirical distribution of the sample median for both cases,
#'   even and odd \eqn{n}, as \code{CDF (x) = pbeta(ECDF(x),(n+1)/2,(n+1)/2)}.
#'
#'   The characteristic function of the empirical distribution of the sample
#'   median is derived as the characteristic function of a weighted discrete
#'   mixture distribution whose support is concentrated on the observed
#'   unique x-values and whose weights are equal to the empirical
#'   probability mass function (PMF) derived from \code{CDF (x) =
#'   pbeta(ECDF(x),(n+1)/2,(n+1)/2)}. In fact, we define the characteristic
#'   function of the sample median as \code{cf = cfE_DiracMixture(t,Xunique,PMF)}.
#'
#'
#'   The characteristic function of the sample median can be used further,
#'   e.g. to test hypotheses about the equality of the medians (see the
#'   Example section).
#'
#' @family Empirical Probability Distribution
#'
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param data vector of observed data  If empty, default value is \deqn{data=0}.
#' @param counts vector of the same size as data that represents the counts
#'            of the repeated values. This is useful if the data are
#'            discrete with many repeted values or if the data are to be
#'            specified from a histogram. If counts is a scalar integer
#'            value, each value specified by data is repeated equally
#'            counts times. If empty, default value is \code{counts = 1}.
#' @param options structure with further optional specificaions:- \code{options.isOgive = false}.
#'
#' @return Characteristic function \eqn{cf(t)} of the empirical distribution of the sample median evalueted at the specified vector argument t.
#' @return Structure with further useful details.
#'
#'   @references
#'  Divine, G. W., Norton, H. J., Barón, A. E., & Juarez-Colunga, E.
#'  (2018). The Wilcoxon–Mann–Whitney procedure fails as a test of
#'  medians. The American Statistician, 72(3), 278-286.
#'
#' @note ver.:02-Sep-2022 18:29:15  (consistent with Matlab CharFunTool, 8-May-2022 17:20:15).
#'
#' @example R/Examples/example_cfE_SampleMedian.R
#'
#' @export
#'
cfE_SampleMedian <- function(t, data, counts, options) {

## CHECK THE INPUT PARAMETERS
if(missing(data)){
        data<-vector()}
if(missing(counts)) {
        counts<-vector()}
if(missing(options)) {
        options<-list()}

if (length(data)==0) {
        data<-0
}
if (length(counts)==0) {
        counts<-1
}

if (is.null(options$isOgive)) {
        options$isOgive <- FALSE
}

l_max <- max(c(length(data), length(counts)))
if (l_max > 1) {
        if (length(data) == 1) {
                data <- rep(data, l_max)
        }
        if (length(counts) == 1) {
                counts <- rep(counts, l_max)
        }

        if ((any(lengths(list(data,counts)) < l_max))) {
                stop("Input size mismatch.")
        }
}

## EMPIRICAL PMF and CDF BASED ON THE OBSERVED DATA
n  <- sum(counts)
nc <- length(counts)
id<-Idx(data)
data <- sort(data)
counts <- counts[id]

id<-firstOccIdx(data)
X <- unique(data)
nX <- length(X)
id <- c(id,nc+1)

Xmin <- min(X)
Xmax <- max(X)
Xdiff <- min(diff(X))

EPMF <- rep(0,length(X))
ECDF <- rep(0,length(X))
cdfsum <- 0

for (i in 1:nX) {
        EPMF[i]<-sum(counts[id[i]:(id[i+1]-1)])/n
        cdfsum<-cdfsum+EPMF[i]
        ECDF[i]<-cdfsum
}



ECDF <- pmax(0,pmin(1,ECDF))

## EMPIRICAL PMF and CDF of the SAMPLE MEDIAN
medianOrder <- (n+1)/2
CDFmedian   <- pbeta(ECDF,medianOrder,medianOrder)
PMFmedian   <- diff(c(0,CDFmedian))

## CHARACTERISTI FUNCTION OF THE SAMPLE MEDIAN
if (options$isOgive){
        if (length(t)==0){
        cf <- function(t) cfE_EmpiricalOgive(t,X,PMFmedian)
        }
        else{
        cf <- cfE_EmpiricalOgive(t,X,PMFmedian)
        }
}
else{
        if (length(t)==0){
        cf <- function(t) cfE_DiracMixture(t,X,PMFmedian)
        }

        else{
        cf <- cfE_DiracMixture(t,X,PMFmedian)
        }
}

## RESULTS
result<-list(
"description" = 'Empirical characteristic function of the sample median',
"cf"          = cf,
"PMFmedian"   = PMFmedian,
"CDFmedian"   = CDFmedian,
"EPMF"        = EPMF,
"ECDF"        = ECDF,
"Xunique"     = X,
"Xmin"        = Xmin,
"Xmax"        = Xmax,
"Xdiff"       = Xdiff,
"data"        = data,
"counts"      = counts,
"n"           = n,
"options"     = options
)
return(result)
}


firstOccIdx <- function(x) {
        x_uniqe <- sort(unique(x))
        indices <- vector()
        for (i in 1:length(x_uniqe)) {
                idx <- 1
                for (j in 1:length(x)) {
                        if (x[j] == x_uniqe[i]) {
                                idx <- j
                                indices <- c(indices, idx)
                                break
                        }
                }
        }
        return(indices)
}

Idx<- function(x) {
        sorted<-sort(x)
        indices <- vector()

        for (i in 1:length(x)) {
                idx<-1
                for (j in 1:length(x)) {
                        is.index<-FALSE
                        if(x[j]==sorted[i]){
                                idx<-j
                                if (length(indices)==0) {
                                        indices<-c(indices,idx)
                                }
                        }
                        if(length(indices)!=0){
                                for (k in 1:length(indices)) {
                                        if(indices[k]==idx){
                                                is.index<-TRUE
                                        }
                                }}
                        if(length(indices)!=0){
                                if(is.index==FALSE){
                                        indices <- c(indices, idx)
                                        break
                                }}


                }

        }
        return(indices)
}

