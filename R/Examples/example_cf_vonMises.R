## EXAMPLE 1:
# CF of the weighted linear combinantion of the von Mises RVs
mu <- c(0, 0, -pi/2, pi/2, 0)
kappa <- c(1, 2, 3, 4, 5)
coef  <- c(1, 2, 3, 4, 5)/15
t <- seq(from = -20,to = 20,length.out = 201)
plotGraf(function(t)
  cf_vonMises(t, mu, kappa, coef), t, title = 'CF of the weighted linear combinantion of the von Mises RVs')

## EXAMPLE2
# CDR/PDF of the weighted linear combinantion of the von Mises RVs
mu <- c(0, 0, -pi/2, pi/2, 0)
kappa <- c(1, 2, 3, 4, 5)
coef  <- c(1, 2, 3, 4, 5)/15
t <- seq(from = -20, to = 20, length.out = 201)
cf <- function(t) {cf_vonMises(t,mu,kappa,coef)}
x <- seq(from = -pi, to = pi, length.out = 201)
prob  <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- -pi
options$xMax <- pi
result <- cf2DistGP(cf,x,prob,options)
result
angle  <- result$x
radius <- result$pdf;
plotPolar(angle, radius)
# figure; polarplot(angle,radius); ???
# x = gca; ax.ThetaAxisUnits = 'radians'; ???

## EXAMPLE3
# CF of the mixture of the von Mises distribution on (-pi,pi)
mu1 <- 0
kappa1 <- 5
mu2 <- 1
kappa2 <- 15
cf <- function(t) {0.25 * cf_vonMises(t, mu1, kappa1) + 0.75 * cf_vonMises(t, mu2, kappa2)}
options <- list()
options$xMin <- -pi
options$xMax <- pi
result <- cf2DistGP(cf = cf, options = options)
angle  <- result$x
radius <- result$pdf
plotPolar(angle, radius)
#   figure; polarplot(angle,radius);
#   ax = gca; ax.ThetaAxisUnits = 'radians';

## EXAMPLE4
# PDF/CDF of the mixture of the von Mises distribution on (0,2*pi)
mu1 <- 0
kappa1 <- 5
mu2 <- 1
kappa2 <- 15
mu3 <- pi
kappa3 <- 10
cf  <- function(t) { 0.25*cf_vonMises(t,mu1,kappa1) +0.25*cf_vonMises(t,mu2,kappa2) + 0.5*cf_vonMises(t,mu3,kappa3)}
options <- list()
options$isCircular   <- TRUE
options$correctedCDF <- TRUE
options$xMin <- 0
options$xMax <- 2*pi
x <- seq(from = 0,to = 2*pi,length.out = 100)
result <- cf2DistGP(cf =cf, x = x, options = options)
angle  <- result$x
radius <- result$pdf
plotPolar(angle, radius)
#figure; polarplot(angle,radius);
#ax = gca; ax.ThetaAxisUnits = 'radians';

