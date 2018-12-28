#' @title Natural Log of the Gamma function
#'
#' @description
#' \code{GammaLog(z)} computes  natural Log of the Gamma function valid in the entire complex
#' plane. This routine uses an excellent Lanczos series
#' approximation for the complex \eqn{ln(Gamma)} function.
#'
#' @family Utility Function
#'
#' @param z may be complex and of any size.
#'
#' @return Natural Log of the Gamma function.
#'
#' @references
#' C. Lanczos, SIAM JNA  1, 1964. pp. 86-96. \cr
#' Y. Luke, "The Special ... approximations", 1969 pp. 29-31. \cr
#' Y. Luke, "Algorithms ... functions", 1977. \cr
#' J. Spouge,  SIAM JNA 31, 1994. pp. 931. \cr
#' W. Press,  "Numerical Recipes". \cr
#' S. Chang, "Computation of special functions", 1996.
#'
#' @note Ver.: 16-Sep-2018 21:15:04 (consistent with Matlab CharFunTool v1.3.0, 24-Jul-2017 10:06:48).
#'
#' @export
#'
GammaLog <- function(z) {
  #CHECK THE INPUT PARAMETERS
  siz <- dim(z)
  zz  <- z

  #f = 0.*z; reserve space in advance

  p <- which(Re(z) < 0)
  if (length(p) > 0) {
    z[p] <- -z[p]
  }

  #ALGORITHM
  #Lanczos approximation for the complex plane

  g <- 607 / 128
  # best results when 4<=g<=5

  c <-
    c(
      0.99999999999999709182,
      57.156235665862923517,
      -59.597960355475491248,
      14.136097974741747174,-0.49191381609762019978,
      0.33994649984811888699e-4,
      0.46523628927048575665e-4,
      -0.98374475304879564677e-4,
      0.15808870322491248884e-3,
      -0.21026444172410488319e-3,
      0.21743961811521264320e-3,
      -0.16431810653676389022e-3,
      0.84418223983852743293e-4,
      -0.26190838401581408670e-4,
      0.36899182659531622704e-5
    )

  s <- 0

  for (k in length(c):2) {
    s <- s + c[k] / (z + (k - 2))
  }

  #for k = size(c,1):-1:2
  #s = s + c(k)./(z+(k-2));
  #end

  zg   <- z + g - 0.5 #vec(zg)
  s2pi <- 0.9189385332046727417803297

  f <- (s2pi + log(c[1] + s)) - zg + (z - 0.5) * log(zg)

  f[z == 1] <- 0.0
  f[z == 2] <- 0.0
  #f[z==1 || z==2] <- 0.0

  if (length(p) > 0) {
    lpi  <- 1.14472988584940017414342735 + 1i * pi
    f[p] <- lpi - log(zz[p] + 0i) - f[p] - log(sin(pi * zz[p]))
  }

  p <- which(round(zz) == zz & Im(zz) == 0 & Re(zz) <= 0)
  if (length(p) > 0) {
    f[p] <- Inf
  }

  if (all(Im(f) == 0)) {
    f <- Re(f)
  }

  dim(f) <- siz

  return(f)
}
