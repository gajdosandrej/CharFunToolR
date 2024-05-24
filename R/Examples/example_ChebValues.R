## EXAMPLE1 (Values of Sine function evaluated Chebyshev points on (-pi,pi))
   n      <- 2^5+1
   domain <- c(-pi,pi)
   x      <- ChebPoints(n,domain)
   f      <-list( sin(x[[1]]))
   coeffs <- ChebCoefficients(f)
   V      <- ChebValues(coeffs)
   print(list(x[[1]], coeffs, f, V))

## EXAMPLE2 (Chebyshev values of the Sine and the Cosine on (-pi,pi))
   n      <- 2^5+1
   domain <- c(-pi,pi)
   x      <- ChebPoints(n,domain)
   f      <-list( sin(x[[1]]), cos(x[[1]]))
   coeffs <- ChebCoefficients(f)
   V      <- ChebValues(coeffs)
   print(list(x[[1]], coeffs, f, V))
