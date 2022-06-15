## EXAMPLE 1
## DISTRIBUTION OF A LINEAR COMBINATION OF THE INDEPENDENT RVs
## (Normal, Student's t, Rectangular, Triangular & Arcsine distribution)
## Y = X_{N} + X_{t} + 5*X_{R} + X_{T} + 10*X_{U}
## CFs: Normal, Student's t, Rectangular, Triangular, and Arcsine
cf_N  <- function(t) exp(-t^2/2)
cf_t <- function(t, nu) {pmin(1, Bessel::BesselK(abs(t) * sqrt(nu), nu / 2, TRUE) *
                                      exp(-abs(t) * sqrt(nu)) *
                                      (sqrt(nu) * abs(t))^(nu / 2) / 2^(nu / 2 - 1)/ gamma(nu / 2))
                        }
#cf_t2 <- function(t) {min(1, Bessel::BesselK(abs(t) * sqrt(1), 1 / 2, TRUE) *
                                                                 #exp(-abs(t) * sqrt(1)) *
                                                                   #   (sqrt(1) * abs(t))^(1 / 2) / 2^(1 / 2 - 1)/ gamma(1 / 2))
#}
cf_R <- function(t) pmin(1, sin(t) / t)
cf_T <- function(t) pmin(1, (2 - 2 * cos(t)) / t^2)
cf_U <- function(t) Bessel::BesselJ(t, 0)
## Characteristic function of the linear combination Y
c <- c(1, 1, 5, 1, 10)
nu <- 1
cf_Y <- function(t) {cf_N(c[1] * t) * cf_t(c[2] * t, nu) * cf_R(c[3] * t) *
                cf_T(c[4] * t) * cf_U(c[5] * t)}
options <- list()
options$N <- 2^10
options$xMin <- -50
options$xMax <- 50
result <- cf2DistFFT(cf = cf_Y, options = options)
# title('CDF of Y = X_{N} + X_{t} + 5*X_{R} + X_{T} + 10*X_{U}')
# problems with BesselK() function...

## EXAMPLE 2
## DISTRIBUTION OF A LINEAR COMBINATION OF THE INDEPENDENT CHI2 RVs
## (Chi-squared RVs with 1 and 10 degrees of freedom)
## Y = 10*X_{Chi2_1} + X_{Chi2_10}
## Characteristic functions of X_{Chi2_1} and X_{Chi2_10}
df1 <- 1
df2 <- 10
cfChi2_1 <- function(t) (1 - 2i * t)^(-df1 / 2)
cfChi2_10 <- function(t) (1 - 2i * t)^(-df2 / 2)
cf_Y <- function(t) cfChi2_1(10 * t) * cfChi2_10(t)
options <- list()
options$xMin <- 0
result <- cf2DistFFT(cf = cf_Y, options = options)
# title('CDF of Y = 10*X_{\chi^2_{1}} + X_{\chi^2_{10}}')

## EXAMPLE3 (PDF/CDF of the compound Poisson-Exponential distribution)
lambda1 <- 10
lambda2 <- 5
cfX <- function(t) cfX_Exponential(t, lambda2)
cf <- function(t) cfN_Poisson(t, lambda1, cfX)
x <- seq(from = 0, to = 8, length.out = 101)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- 1
result <- cf2DistFFT(cf, x, prob, options)
