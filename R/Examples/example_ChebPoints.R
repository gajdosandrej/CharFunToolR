## EXAMPLE 1
# Barycentric interpolant of the Sine function on (-pi,pi)
x <- ChebPoints(32, c(-pi, pi))
f <- sin(x)
xNew <- seq(from = -pi,
            to = pi,
            length.out = 201)
fNew <- interpBarycentric(x, f, xNew)[[2]]
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
lines(xNew, fNew, col = "blue",  lwd = 2)
print(list(
        "xNew" = xNew,
        "fNew" = fNew,
        "sin(xNew)" = sin(xNew)
))
