## FUNCTION

FindRoots  <- function(fun,A,B,n,isplot,tol,maxiter){

## CHECK THE INPUT PARAMETERS

if (missing(maxiter)) {
        maxiter <- vector()
}

if (missing(tol)){
        tol <- vector()
}
if (missing(isplot)){
        isplot <- vector()
}
if (missing(n)){
        n <- vector()
        }


if (missing(B)){
        B <- vector()
        }
if (missing(A)){
        A <- vector()
}

## SET THE DEFAULT VALUES (input parameters)

if (length(maxiter)==0){
maxiter <- 100
}

if (length(tol)==0){
tol <- 1e-16
}

if (length(isplot)==0){
isplot <- TRUE
}

if (length(n)==0){
n <- 2^5
}

if (length(B)==0){
B <- 1
}

if (length(A)==0){
A = -1
}

interval <- matrix(c(A,B))

if (abs(B-A) > 3*pi){
divisionrule <- c(-0.5,0,0.5)
interval     <-  GetSubs(interval,divisionrule)
}

else{
        divisionrule <- 0
}

## ALGORITHM

cgl <- setupChebyshev(n)$cgl
M   <- setupChebyshev(n)$M

roots   <- vector()
warning <- 0
iter    <- 0

while (iter < maxiter){
        iter        <- iter + 1
        result      <- rootFinder(cgl,M,fun,interval,n,tol)
        r           <- result$roots
        interval    <- result$intervals
        err         <- result$err
        isWarning   <- result$isWarning
        roots       <- c(roots,r)

        if (length(interval)==0){
                break
        }
        else{
        interval <-  GetSubs(interval,divisionrule)

        }
        warning  <- warning + isWarning
}

if (iter == maxiter){
warning <- warning + 1
}

roots  <- unique(roots)
result <- list("roots"=roots, "warning"=warning,"err"=err)

#if (isplot){
#N <- 1000
#t <- seq(A,B,length.out=N)
#plot(t,fun(t))

#plot(roots,fun(roots))

#}
return(result)
}

## Function setupChebyshev
setupChebyshev <- function (n){
# setupChebyshev evaluates the Chebyshev-Gauss-Legendre points and the
# Chebyshev Transformation Matrix
#
# Adapted from the MATLAB code published in Day & Romero (2005): Roots of
# polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
# Numerical Analysis,  43, 1969 - 1987.

# Chebyshev-Gauss-Legendre points
ind <- c(0:n)
cgl <- cos(ind * pi/n)

# Chebyshev Transformation Matrix
M        <- cos(ind%*% t(ind) * pi/n)
M[1, ]   <- M[1, ]/2
M[ ,1]   <- M[ ,1]/2
M[n+1, ] <- M[n+1, ]/2
M[ ,n+1] <- M[ ,n+1]/2
M        <- M * (2/n)

result  <- list(
        "cgl" = cgl,
        "M"   = M)

return(result)
}


## Function evalCheb
evalCheb <- function (n, z) {
#EVALCHEB evaluates the required Vandermonde matrix
#
# Input: vector of points, z, and the polynomial degree n.
# Output: Vandermonde matrix, m by degree_max + 1, V(j+1,k)=T_j(z_k)
#
# Adapted from the MATLAB code published in Day & Romero (2005): Roots of
# polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
# Numerical Analysis,  43, 1969 - 1987.

        m <- length(z)
        V <- matrix()
        if (m*n >= 0){
                V <- cbind(rep(1, m))
                if (n >= 1){
                        V <- cbind(V,z)
                        if (n >= 2){
                                index <- pracma::finds(log(abs(z)) >= 100/n )
                                si    <- length(index)
                                if (si > 0){
                                        z[index] <- rep(1, si)*exp(100/n)
                                }
                                for (i in 2:n){
                                        V <- cbind(V,V[ ,i]*(2*z) - V[ ,i-1])
                                }
                        }
                }
        } else {

                V <- vector()
        }


        return(V)
}

rootFinder <- function(cgl,M,fun,intervals,n,tol){
#ROOTFINDER estimates the roots of fun over intervals by using nth order
# Chebyshev polynomial approximation of fun
#
# Adapted from the MATLAB code published in Day & Romero (2005): Roots of
# polynomials expressed in terms of orthogonal polynomials. SIAM Journal on
# Numerical Analysis,  43, 1969 - 1987.

# Ver.: Sat Feb 11 15:14:05 2023

roots <- vector()

nint <- dim(intervals)[2]
err  <- matrix(FALSE,nrow = nint, ncol = 1)
isWarning <- 0
for (i in 1:nint){
        range <- (intervals[2,i] - intervals[1,i])/2
        shift <- (intervals[1,i] + intervals[2,i])/2
        x <- cgl*range + shift
        f <- fun(x)
        ExpansionCoeff <- t(f) %*% M
        if (abs(ExpansionCoeff[n+1]) < tol){
        isWarning <- 1
        #warning('The leading expansion coefficient vanishes')

        } else {
                ExpansionCoeff <- ExpansionCoeff/(-2*ExpansionCoeff[n+1])
                H      <- mrbsizeR::tridiag(0.0,rep(1, n-1)/2, 0.0) + mrbsizeR::tridiag(0,0,rep(1, n-1)/2)
                H[1,2] <- 1
                C      <- H
                C[n, ] <- C[n, ] + ExpansionCoeff[1:n]
                Eigenvalues <- pracma::eig(C)
                Vandermonde <- evalCheb(n,Eigenvalues)
                Vcolsums    <- colSums(abs(t(Vandermonde)))
                tube_index  <- pracma::finds((abs(Im(Eigenvalues))<2)
                       && (abs(Re(Eigenvalues))< 2))
                Solutions   <- Eigenvalues[tube_index]
                Vcolsums    <- Vcolsums[tube_index]
                cond_max    <- min(2^(n/2), 10^6)
                condEigs_index <-  Vcolsums < cond_max
                Solutions   <- sort(Solutions[condEigs_index])
                r <- range * Solutions + shift
                if (!is.complex(r)) {
                        r     <- r[r>=intervals[1,i] && r<=intervals[2,i]]
                        roots <- rbind(roots,r)
                } else {
                        err[i] <- TRUE
                }
          }
        }
intervals <- matrix(intervals[ ,err],nrow=length(intervals[,1]))

result <- list(
        "roots"=roots,
        "intervals"=intervals,
        "err"=err,
        "isWarning"=isWarning)

return(result)
}

## FUNCTION GETSUBS (Sub-division of the integration intervals)
GetSubs <-function (SUBS, XK)  {
#GETSUBS Sub-division of the intervals for adaptive root finding
#
# Ver.: Sat Feb 11 15:57:13 2023

NK <- length(XK)
SUBIND <- c(1:NK)
M <- 0.5*(SUBS[2, ]-SUBS[1, ])
C <- 0.5*(SUBS[2, ]+SUBS[1, ])
Z<-matrix()

Z <- rbind(XK*M + rep(1, NK)*C)

L <- c(SUBS[1, ], Z[SUBIND, ])
U <- c(Z[SUBIND, ], SUBS[2, ])
SUBS <- matrix(c(matrix(L, nrow = 1), matrix(U, nrow = 1)),nrow = 2,byrow = TRUE)
return(SUBS)
}

