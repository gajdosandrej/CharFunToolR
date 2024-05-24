## EXAMPLE 1
# Barycentric interpolant of the Sine function on (-pi,pi)
domain<-c(-pi,pi)
x <- ChebPoints(32, c(-pi, pi))
f <- list(sin(x[[1]]))
xNew <- seq(from = -pi,
            to = pi,
            length.out = 201)
fNew <- interpBarycentric(x[[1]], f[[1]], xNew)[[2]]
fCheb<-InterpChebValues(f,xNew,domain)
fTrue<-sin(xNew)
plot(
        x[[1]],
        f[[1]],
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
        "fTrue"=fTrue,
        "fCheb"=fCheb,
        "xNew" = xNew,
        "fNew" = fNew,
        "sin(xNew)" = sin(xNew)
))

## EXAMPLE 2
#Integral of the Sine function on the interval (0,pi)
domain<-c(0,pi)
x<-ChebPoints(32,domain)
f<-list(sin(x[[1]]))
Itrue<-2
Icalc<-t(x[[2]])%*%f[[1]]
print(list(
        "Itrue"=Itrue,
        "Icalc"=Icalc
))

