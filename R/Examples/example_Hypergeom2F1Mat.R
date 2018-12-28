## EXAMPLE
a <- 3
b <- 2.5
c <- 1.5
# X <- c(1, 2, 3) / 5
X <- t(c(1, 2, 3) / 5)
MAX <- 50
f <- Hypergeom2F1Mat(a, b, c, X, MAX)

# a <- 3
# b <- c(1,2,3,4,5)
# c <- c(5,4,3,2,1)
# X <-  t(c(1, 2, 3)) / 5
# MAX <- 10
# f <- Hypergeom2F1Mat(a, b, c, X, MAX)

