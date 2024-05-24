## EXAMPLE 1
# CF of the symmetric Wigner distribution on (-1,1)
t <- seq(from = -50,
         to = 50,
         length.out =501)
plotReIm(function(t)
        cfS_Wigner(t),
        t,
        title = "CF of the Wigner distribution on (-1,1)")



##EXAMPLE2
# PDF/CDF of Wigner distribution on (-1,1)
cf <- function(t)
        cfS_Wigner(t)
x <- seq(-1,1,length.out = 501)
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$xMin <- -1
options$xMax <- 1
result <- cf2DistGP(cf, x, prob, options)

