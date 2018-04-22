## EXAMPLE 1
# CF of the Rectangular distribution on (-1,1)
t <- seq(from = -50,to = 50,length.out = 201)
plotGraf(function(t)
  cf_RectangularSymmetric(t), t, title = "CF of the Rectangular distribution on (-1,1)")

## EXAMPLE 2
# CF of a weighted linear combination of independent Rectangular RVs
t <- seq(from = -10,to = 10,length.out = 201)
coef <- c(1, 2, 3, 4, 5)/15
plotGraf(function(t)
  cf_RectangularSymmetric(t,coef), t, title = "CF of a weighted linear combination of Rectangular RVs")

## EXAMPLE 3
# PDF/CDF of a weighted linear combination of independent Rectangular RVs
coef <- c(1, 2, 3, 4, 5)/15
cf <- function(t) cf_RectangularSymmetric(t,coef)
x <- seq(from = -1,to = 1,length.out = 201)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$N <- 2^12
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf,x,prob,options)
