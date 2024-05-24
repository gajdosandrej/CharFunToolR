#' @title Empirical characteristic function of the sample mean based on the observed random sample
#'
#' @description
#' The empirical distribution of the sample mean is estimated by the
#' bootstrap distribution of the sample mean (i.e. the bootstrap mean
#' distribution).
#'
#'   The characteristic function of the empirical distribution of the sample
#'   mean is derived as the characteristic function of an equally weighted
#'   linear combination (arithmetic mean) of n independent random variables
#'   each having distribution defined by the observed empirical cumulative
#'   distribution function \eqn{(ECDF)}, defined from the observed data. This can
#'   be expressed as a weighted mixture of the DIRAC distributions whose
#'   support is concentrated on the observed unique x-values and the weights
#'   are equal to the empirical probability mass function \eqn{(EPMF)} derived
#'   from the \eqn{ECDF}.
#'
#' In fact, we define the characteristic function of the sample mean as \code{cf =cfE_DiracMixture(t/n,X,EPMF)^n)}, where
#' \eqn{X} is a vector of the observed unique x-values and \eqn{EPMF} is the vector of the weights.
#'
#' The characteristic function of the sample mean can be used further,
#' e.g. to test hypotheses about the equality of the means (see the
#' Example section).
#'
#'
#' @family Empirical Probability Distribution
#'

#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param data vector of original data. If empty, default value is \code{data = 0}.
#' @param counts vector of the same size as data that represents the counts
#'            of the repeated values. This is useful if the data are
#'            discrete with many repeted values or if the data are to be
#'            specified from a histogram. If counts is a scalar integer
#'            value, each value specified by data is repeated equally
#'            counts times. If empty, default value is \code{counts = 1}.
#' @param options structure with further optional specificaions:- \code{options.isOgive = false}.
#'
#' @return Characteristic function \eqn{cf(t)} of the empirical distribution of the
#' sample mean evalueted at the specified vector argument t.
#'
#' @return  structure with further useful details.
#'
#' @note Ver.: 02-Sep-2022 17:53:40 (consistent with Matlab CharFunTool, 08-May-2022 17:20:15).
#'
#' @example R/Examples/example_cfE_SampleMean.R
#'
#' @export
#'


cfE_SampleMean <- function(t, data, counts, options) {
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


        ##EMPIRICAL PMF and CDF BASED ON THE OBSERVED DATA

        n  <- sum(counts)
        nc <- length(counts)
        id<-Idx(data)
        data <- sort(data)
        counts <- counts[id]

        id<-firstOccIdx(data)
        X <- unique(data)
        nX <- length(X)
        id <-c(id,nc+1)

        Xmin <- min(X)
        Xmax <- max(X)
        Xdiff <- min(diff(X))



        Xcounts <- rep(0,length(X))
        EPMF <- rep(0,length(X))
        ECDF <- rep(0,length(X))
        cdfsum <- 0

        for (i in 1:nX) {
        Xcounts[i]<-sum(counts[id[i]:(id[i+1]-1)])
        EPMF[i]<-Xcounts[i]/n
        cdfsum<-cdfsum+EPMF[i]
        ECDF[i]<-cdfsum
        }

        ECDF <- pmax(0,pmin(1,ECDF))

        ## CHARACTERISTI FUNCTION OF THE SAMPLE MEAN

if (options$isOgive){

        if(length(t)==0 ){
                cf <- function(t) cfE_EmpiricalOgive(t/n,X,EPMF)^n

        }

         else{
                cf <- cfE_EmpiricalOgive(t/n,X,EPMF)^n
         }
}
else{

         if (length(t)==0){
        cf <- function(t) cfE_DiracMixture(t/n,X,EPMF)^n
        }
        else{
        cf <- cfE_DiracMixture(t/n,X,EPMF)^n
        }

}

        ## RESULTS
        result<-list(
        "description" = 'Empirical characteristic function of the sample mean',
        "cf"          = cf,
        "EPMF"        = EPMF,
        "ECDF"        = ECDF,
        "X"           = X,
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

