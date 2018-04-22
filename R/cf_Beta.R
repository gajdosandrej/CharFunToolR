#' @title Characteristic function of a linear combination of independent BETA random variables
#'
#' @description
#' \code{cf_Beta(t, alpha, beta, coef, niid)} evaluates the characteristic function of a linear combination
#' (resp. convolution) of independent BETA random variables defined on the interval \eqn{(0,1)}.
#'
#' That is, \code{cf_Beta} evaluates the characteristic function \eqn{cf(t)}
#' of  \eqn{Y = sum_{i=1}^N coef_i * X_i}, where \eqn{X_i ~ Beta(\alpha_i,\beta_i)}
#' are independent RVs, with the shape parameters \eqn{\alpha_i > 0} and \eqn{\beta_i >0},
#' and with the \eqn{mean = \alpha_i / (\alpha_i + \beta)_i} and the
#' \eqn{variance = (\alpha_i*\beta_i) / ((\alpha_i+\beta_i)^2*(\alpha_i+\beta_i+1))}, for \eqn{i = 1,...,N}.
#'
#' The characteristic function of \eqn{X ~ Beta(\alpha,\beta)} is
#' \deqn{cf(t) = cf_Beta(t,\alpha,\beta) = 1F1(\alpha; \alpha + \beta; i*t),}
#' where \eqn{1F1(.;.;.)} is the Confluent hypergeometric function. Hence,
#' the characteristic function of \eqn{Y  = coef(1)*X_1 + ... + coef(N)*X_N}
#' is \eqn{cf(t) =  cf_X_1(coef(1)*t) * ... * cf_X_N(coef(N)*t)},
#' where \eqn{X_i ~ Beta(\alpha(i),\beta(i))} with \eqn{cf_X_i(t)}.
#'
#' @family Continuous Probability distribution
#'
#' @seealso For more details see WIKIPEDIA:
#' \url{https://en.wikipedia.org/wiki/Beta_distribution}.
#'

#' @param t vector or array of real values, where the CF is evaluated.
#' @param alpha vector of the 'shape' parameters \code{alpha > 0}. If empty, default value is \code{alpha = 1}.
#' @param beta vector of the 'shape' parameters \code{beta > 0}. If empty, default value is \code{beta = 1}.
#' @param coef vector of the coefficients of the linear combination of the Beta distributed random variables.
#' If coef is scalar, it is assumed that all coefficients are equal. If empty, default value is \code{coef = 1}.
#' @param niid scalar convolution coeficient \code{niid}, such that \eqn{Z = Y + ... + Y}
#' is sum of \eqn{niid} iid random variables \eqn{Y}, where each \eqn{Y = sum_{i=1}^N coef(i) * log(X_i)}
#' is independently and identically distributed random variable. If empty, default value is \code{niid = 1}.
#'
#' @return Characteristic function \eqn{cf(t)} of a linear combination of independent BETA random variables.
#'
#' @example R/Examples/example_cf_Beta.R
#'
#' @export
#'
cf_Beta <- function(t, alpha, beta, coef, niid) {
  ## CHECK THE INPUT PARAMETERS
  if(missing(alpha)) {
    alpha <- vector()
  }
  if(missing(beta)) {
    beta <- vector()
  }
  if(missing(coef)) {
    coef <- vector()
  }
  if(missing(niid)) {
    niid <- vector()
  }

  if(length(beta) == 0 && length(alpha) > 0) {
    beta <- 1
  } else if(length(beta) == 0 && length(coef) > 0) {
    beta <- 1
  } else if(any(beta == 0) || length(beta) == 0){
    beta <- 1
  }

  if(length(alpha) == 0 && length(beta) > 0) {
    alpha <- 1
  } else if(length(alpha) == 0 && length(coef) > 0) {
    alpha <- 1
  }

  if(length(coef) == 0 && length(alpha) > 0) {
    coef <- 1
  } else if(length(coef) == 0 && length(beta) > 0) {
    coef <- 1
  }

  if(length(niid) == 0) {
    niid <- 1
  }

  ## Equal size of the parameters
  if(length(coef) > 0 && length(alpha) == 1 && length(beta) == 1 && length(niid) == 0) {
    coef <- sort(coef)
    coef_orig_sort <- coef
    m <- length(coef)
    coef <- unique(coef)
    idx <- firstOccIdx(coef_orig_sort)
    alpha <- alpha * diff(c(idx,m+1))
  }

  l_max <- max(c(length(coef), length(alpha), length(beta)))
  if (l_max > 1) {
    if (length(alpha) == 1) {
      alpha <- rep(alpha, l_max)
    }
    if (length(coef) == 1) {
      coef <- rep(coef, l_max)
    }
    if (length(beta) == 1) {
      beta <- rep(beta, l_max)
    }
    if ((any(lengths(list(coef, alpha, beta)) < l_max))) {
      stop("Input size mismatch.")
    }
  }

  ## Characteristic function
  szc <- length(coef)
  szt <- dim(t)
  t <- c(t)

  cf <- 1
  for(i in 1:szc) {
    cf <- cf * hypergeom1F1(alpha[i],alpha[i]+beta[i],1i*coef[i]*t)$f
  }

  dim(cf) <- szt
  cf[t == 0] <- 1

  if (length(niid) > 0) {
    if (length(niid) == 1) {
      cf <- cf ^ niid
    } else {
      stop("niid should be a scalar (positive integer) value.")
    }
  }

  return(cf)

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
