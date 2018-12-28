## EXAMPLE 1
n <- 0:5
x <- ChebPoints(32, c(-1, 1))
pval <- ChebPoly(n, x)
matplot(x, t(pval), main = 'Chebyshev polynomials T_n(x) evaluated at Chebyshev points', type = "b")
