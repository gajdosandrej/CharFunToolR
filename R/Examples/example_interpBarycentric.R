## EXAMPLE 1
# Barycentric interpolant of the Sine function on (-pi,pi)
x <- ChebPoints(32, c(-pi, pi))[[1]]
f <- sin(x)
xNew <- seq(from = -pi,
            to = pi,
            length.out = 201)
fNew <- interpBarycentric(x, f, xNew)
plot(
        x,
        f,
        type = "p",
        xlab = "",
        ylab = "",
        main = "",
        pch = 20,
        col = "red",
        cex = 1.5
)
lines(xNew, fNew[[2]], col = "blue",  lwd = 2)
print(list(
        "xNew" = xNew,
        "fNew" = fNew[[2]],
        "sin(xNew)" = sin(xNew)
))

## EXAMPLE 2
# Barycentric interpolant of the Normal CDF
x <- ChebPoints(21, c(-8, 8))
f <- pnorm(x[[1]])
xNew <- seq(from = -5,
            to = 5,
            length.out = 201)
fNew <- interpBarycentric(x[[1]], f, xNew)[[2]]
plot(
        xNew,
        fNew,
        type = "p",
        xlim = c(-8 , 8),
        xlab = "",
        ylab = "",
        main = "",
        pch = 20,
        col = "red",
        cex = 1
)
lines(x[[1]], f, col = "blue",  lwd = 2)
print(list(
        "xNew" = xNew,
        "fNew" = fNew,
        "pnorm(xNew)" = pnorm(xNew)
))

## EXAMPLE 3
# Barycentric interpolant of the Normal quantile function
x <- ChebPoints(2 ^ 7, c(-8, 8))
cdf <- pnorm(x[[1]])
prob <- seq(from = 1e-4,
            to = 1 - 1e-4,
            length.out = 101)
qf <- interpBarycentric(cdf, x[[1]], prob)[[2]]
plot(
        cdf,
        x[[1]],
        type = "p",
        xlim = c(0, 1),
        xlab = "",
        ylab = "",
        main = "",
        pch = 20,
        col = "red",
        cex = 1
)
lines(prob, qf, col = "blue",  lwd = 2)
print(list(prob, qf, qnorm(prob)))

## EXAMPLE 4
# Barycentric interpolant of the Compound distribution
cfX <- function(t)
        cfX_Exponential(t, 5)
cf <- function(t)
        cfN_Poisson(t, 10, cfX)
x <- ChebPoints(2 ^ 9, c(0, 10))
prob <- c(0.9, 0.95, 0.99)
options <- list()
options$isCompound <- TRUE
result <- cf2DistGP(cf, x[[1]], prob, options)
CDF <- function(x)
        interpBarycentric(result$x, result$cdf, x)[[2]]
PDF <- function(x)
        interpBarycentric(result$x, result$pdf, x)[[2]]
QF <- function(p)
        interpBarycentric(result$cdf, result$x, p)[[2]]
prob <- seq(from = 1e-4,
            to = 1 - 1e-4,
            length.out = 201)
plot(
        prob,
        QF(prob),
        type = "p",
        xlim = c(0, 1),
        xlab = "",
        ylab = "",
        main = "",
        pch = 20,
        col = "red",
        cex = 1
)
lines(result$cdf, result$x, col = "blue",  lwd = 2)
