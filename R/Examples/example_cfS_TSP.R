## EXAMPLE 1
# CF of the symmetric TSP distribution with theta =3/2 on (-1,1)
theta <- 3/2
t <- seq(from = -50,
         to = 50,
         length.out =501)
plotReIm(function(t)
        cfS_TSP(t, theta),
        t,
        title = "CF of the symmetric TSP distribution  on (-1,1)")



##EXAMPLE2
# PDF/CDF of the symmetric TSP distribution on (-1,1)
thet <- 3/2
cf <- function(t)
        cfS_TSP(t, theta)
x <- seq(-1,1,length.out = 101)
xRange <- 2
options <- list()
options.N <-2^8
options.dt <- 2*pi/xRange
result <- cf2DistGP(cf, x, c(), options)
