#' @title Characteristic function of the Ogive empirical distribution based on the observed histogram
#'
#' @description
#' Characteristic function of the OGIVE EMPIRICAL distribution based on
#' the observed histogram given by data and counts (the vector lengths of
#' data and counts are equal) or by edges and counts (the length of the
#' edges is equal to the length of counts + 1).
#'
#' That is, \eqn{cf(t)} is given as a weighted mixture of the RECTANGULAR
#' (continuous uniform) RVs,
#' \code{cf(t) = cfE_EmpiricalOgive(t,data,counts) = sum_{j=1}^n freq_j * cf_Rectangular(t,edges_{j},edges_{j+1})},
#' where \code{cf_Rectangular(t,edges_{j},edges_{j+1})} is CF of the RECTANGULAR
#' distribution on the interval \eqn{[edges_{j},edges_{j+1}]}. Here, the
#' frequencies are derived from the counts, \eqn{FREQ = COUNTS/SUM(COUNTS)}, and
#' the vector EDGES is derived automatically from the DATA (or is
#' specified by the user as an input vector EDGES of size n+1, instead of
#' DATA). In particular, for \deqn{j = 2,...,n} we set \eqn{EDGES(j) = (DATA(j-1) +
#' DATA(j))/2} and further,  \eqn{EDGES(1) = DATA(1) - (DATA(2)-DATA(1))/2} and
#' \code{EDGES(n+1) = DATA(n) + (DATA(n)-DATA(n-1))/2}.
#'
#' In particular, \eqn{[COUNTS,EDGES] = HISTCOUNTS(DATA)} partitions the values
#' in DATA into bins, and returns the count in each bin, as well as the
#' bin edges. HISTCOUNTS uses an automatic binning algorithm that returns
#' bins with a uniform width, chosen to cover the range of elements in
#' DATA and reveal the underlying shape of the distribution. \eqn{COUNTS(k)}
#' will count the value \eqn{DATA(i)} if \deqn{EDGES(k) <= DATA(i) < EDGES(k+1)}. The
#' last bin will also include the right edge such that COUNTS(end) will
#' count \eqn{DATA(i)} if \deqn{EDGES(end-1) <= DATA(i) <= EDGES(end)}.
#'
#'
#' @family Empirical Probability Distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Empirical_distribution_function}.
#' \url{https://en.wikipedia.org/wiki/Ogive_(statistics)}
#'
#' @param t vector or array of real values, where the CF is evaluated.
#' @param data vector of data (if the size is equal to the vector counts) or a vector of edges (if)the size is \eqn{n+1}, where \eqn{n} is the length of counts). If empty, default value is \deqn{data=1}.
#' @param counts vector of counts (integers) related to the values given by data. If empty, default value is \deqn{couns=1}.
#'
#' @return Characteristic function \eqn{cf(t)} of the Ogive empirical distribution, based on the observed data.
#'
#' @references
#' WITKOVSKY V., WIMMER G., DUBY T. (2017). Computing the aggregate
#' loss distribution based on numerical inversion of the compound empirical
#' characteristic function of frequency and severity. arXiv preprint arXiv:1701.08299.
#'
#' @note ver.:18-aug-2022 17:50:50  (consistent with Matlab CharFunTool, 8-May-2022 15:45:27).
#'
#' @example R/Examples/example_cfE_EmpiricalOgive.R
#'
#' @export

cfE_EmpiricalOgive <- function(t, data, counts) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(data)){
                data <- vector()}
        if(missing(counts)) {
                counts <- vector()}

        if (length(data)==0) {
                data <- 1
        }

        if (length(counts)==0) {
                counts <- 1
        }

        if(length(counts)==1){

                l_max <- max(c(length(data), length(counts)))
                if (l_max > 1) {
                        if (length(data()) == 1) {
                                data() <- rep(data, l_max)
                        }
                        if (length(counts) == 1) {
                                counts <- rep(counts, l_max)
                        }
                        if (any(lengths(list(data, counts)) < l_max)) {
                                stop("Input size mismatch.")
                        }
                }
        }



        data<-Conj(data)
        counts<-Conj(counts)

        n<-length(counts)
        freq<-counts/sum(counts)
   if(length(data)==n){

        if (length(data)==1){
           edges<-rep(0,2)
           dspan<-1
           edges[1]<-data[1]-dspan
           edges[2]<-data[1]+dspan
        }
     else if (length(data)==2){
          edges<-rep(0,3)
          dspan<-(data[2]-data[1])/2
          edges[1]<-data[1]-dspan
          edges[2]<-data[1]+dspan
          edges[3]<-data[2]+dspan

        }
        else{
        edges<-rep(0,n+1)
        edges[1]<-data[1]-(data[2]-data[1])/2
        edges[2:n]<-(data[1:(n-1)]+data[2:n])/2
        edges[n+1]<-data[n]+(data[n]-data[n-1])/2

}
   }

 else if (length(data)==n+1){
                edges<-data

        }

        ## Characteristic function
        szt <- dim(t)
        t <- t(t)
        t <- t(t)
        cf <- 0

        aux <- exp(1i * t %*% edges)
        aux <- (aux[ ,2:(n+1)]- aux[ ,1:n]) / (1i * t %*% (edges[2:(n+1)]-edges[1:n]))

                cf <- cf + rowSums(freq * aux)

        cf[t==0] <- 1
        dim(cf) <- szt

        return(cf)
}
