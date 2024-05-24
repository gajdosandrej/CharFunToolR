## EXAMPLE 1
# CF of the distribution of ChiSquare RV with DF = 1, NCP = 1
df <- 1
ncp <- 1
coef <- 1
t <- seq(from = -100,
         to = 100,
         length.out = 501)
plotReIm(
        function(t)
                cf_ChiSquare(t, df, ncp, coef),
        t,
        title = expression('Characteristic function of the' ~ chi ^ 2 ~ 'RV with DF=1 and NCP=1')
)

## EXAMPLE 2
# CF of a linear combination of K=100 independent ChiSquare RVs
coef <- 1 / (((1:50) - 0.5) * pi) ^ 2
plot(
        1:50,
        coef,
        xlab = "",
        ylab = "",
        type = "p",
        pch = 20,
        col = "blue",
        cex = 1,
        main = expression('Coefficients of the linear combination of' ~ chi ^ 2 ~ 'RVs with DF=1')
)
lines(1:50, coef, col = "blue")
df <- 1
t <- seq(from = -100,
         to = 100,
         length.out = 501)
plotReIm(
        function(t)
                cf_ChiSquare(t = t, df = df, coef = coef),
        t,
        title = expression('CF of the linear combination of' ~ chi ^ 2 ~ 'RVs with DF=1')
)

## EXAMPLE 3
# PDF/CDF from the CF by cf2DistGP
df <- 1
coef <-  1 / (((1:50) - 0.5) * pi) ^ 2
cf <- function(t)
        cf_ChiSquare(t = t, df = df, coef = coef)
options <- list()
options$N <- 2 ^ 12
options$xMin <- 0
x <- seq(from = 0,
         to = 3.5,
         length.out = 501)
prob <- c(0.9, 0.95, 0.975, 0.99)
result <- cf2DistGP(cf, x, prob, options)
result
