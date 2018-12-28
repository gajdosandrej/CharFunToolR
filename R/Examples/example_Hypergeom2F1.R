## EXAMPLE 1
a <- 3
b <- 2.5
c <- 1.5
z <- 1i * seq(0, 1, by = 0.05)
f <- Hypergeom2F1(a, b, c, z)

# EXAMPLE 2
t <- 1i * seq(-5, 5, length.out = 11)
a <- 3 * t
b <- 2.5 * t
c <- 1.5 * t
z <- 0.75
f <- Hypergeom2F1(a,b,c,z)
