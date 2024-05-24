#' @title
#' Hypergeometric function of a matrix argument
#'
#' @description
#' \code{HypergeompFqMat(a, b, x, y, alpha, MAX, lam)} computes the truncated hypergeometric
#' function \eqn{pFq^alpha(a;b;x;y)} with parameters \code{a} and \code{b} and of two matrix arguments, \code{x} and \code{y}.
#'
#' @details
#' Here, hypergeometric function \eqn{pFq^alpha(a;b;x;y)} is defined as
#' \deqn{pFq^alpha(a;b;x;y) = pFq^alpha(a1,...,ap;b1,...,bp;x;y) =
#' sum_k sum_kappa [((a1)_kappa^(alpha) ... (a1)_kappa^(alpha)) /
#' (k!(b1)_kappa^(alpha) ... (b1)_kappa^(alpha)] C_kappa^alpha(X),}
#' where \eqn{kappa = (kappa1,kappa2,...)} is a partition of \eqn{k} and \eqn{(ak)_kappa^(alpha), (bk)_kappa^(alpha)}
#' denote the generalized Pochhammer symbols, and \eqn{C_kappa^alpha(X)} is the Jack function.
#'
#' For statistical applications considered here we need to use the confluent hypergeometric function
#' \eqn{1F1(a;b;X)} and the Gauss hypergeometric function \eqn{2F1(a,b;c;X)} in matrix argument.
#' These can be expressed as special cases of the generalized hypergeometric function \eqn{pFq^alpha(a;b;X)},
#' with the parameter \eqn{alpha = 2} (the case with zonal polynomials).
#'
#' For more details and definition of the hypergeometric function \eqn{pFq^alpha(a;b;x;y)} with matrix argument see,
#' e.g., Koev and Edelman (2006) or Muirhead (2009).
#'
#' @author
#' Copyright (c) 2004 Plamen Koev.
#' MHG is a collection of MATLAB functions written in C for computing the hypergeometric function
#' of a matrix argument.
#'
#' MHG LICENCE: \cr
#' This program is free software; you can redistribute it and/or modify it
#' under the terms of the GNU General Public License as published by the Free
#' Software Foundation; either version 2 of the License, or (at your option) any later version.
#'
#' This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#' without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#' See the GNU General Public License for more details.
#'
#' You should have received a copy of the GNU General Public License along with this program;
#' if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#'
#' If you use MHG in any program or publication, please acknowledge its author by adding a reference.
#'
#' This is R version of modified (by Viktor Witkovsky) version of the original MATLAB code hg.m by Plamen Koev
#'
#' @references
#' [1] Koev, P. and Edelman, A., 2006. The efficient evaluation of the hypergeometric function
#' of a matrix argument. \emph{Mathematics of Computation}, 75(254), 833-846.
#'
#' [2] Muirhead RJ. Aspects of multivariate statistical theory. John Wiley & Sons; 2009 Sep 25.
#'
#' @param a parameter of the hypergeometric function \eqn{pFq^alpha(a;b;x;y)}.
#' @param b parameter of the hypergeometric function \eqn{pFq^alpha(a;b;x;y)}.
#' @param x matrix argument (given as vector of eigenvalues).
#' @param y optional second matrix argument (vector of eigenvalues).
#' @param alpha parameter of the hypergeometric function \eqn{pFq^alpha(a;b;x;y)}, default value is \code{alpha = 2}.
#' @param MAX maximum number of partitions, \eqn{|kappa|<=MAX}, default value is \code{MAX = 20}.
#' @param lam optional parameter, \eqn{kappa<=lam}.
#'
#' @family Utility Function
#'
#' @return
#' List including the following items:
#'
#' \item{s}{hypergeometric sum, \eqn{pFq^alpha(a;b;x;y)},}
#' \item{ss}{partial sums.}
#'
#' @note Ver.: 10-Dec-2018 17:35:21 (consistent with Matlab CharFunTool v1.3.0, Ver.: 23-Oct-2017 11:47:05).
#'
#' @example R/Examples/example_HypergeompFqMat.R
#'
#' @export
#'
HypergeompFqMat <- function(a, b, x, y, alpha = 2, MAX = 20, lam) {

        ## CHECK THE INPUT PARAMETERS
        if(missing(y)) {
                y <- vector()
        }
        if(missing(lam)) {
                lam <- vector()
        }

        ## CHECK THE COMMON SIZE of the parameters a and b
        # if(is.matrix(a) && is.matrix(b) && nrow(a) != nrow(b)) {
        #         stop('Input size mismatch.')
        # }
        # if(is.vector(a) && is.vector(b) && length(a) != length(b)) {
        #         stop('Input size mismatch.')
        # }

        na <- numeric()
        nb <- numeric()
        if(is.vector(a)) {
                na <- length(a)
        } else if (is.matrix(a)) {
                na <- nrow(a)
        }
        if(is.vector(b)) {
                nb <- length(b)
        } else if (is.matrix(b)) {
                nb <- nrow(b)
        }

        ## CHECK THE COMMON SIZE of the parameters a and b
        if(na != nb) {
                stop('Input size mismatch.')
        }

        ## ALGORITHM
        if(is.vector(x)) {
                n <- length(x)
        } else if(is.matrix(x)) {
                n <- max(dim(x))
        }

        lambda <- floor(MAX/(1:n))

        if(length(lam) > 0) {
                lambda <- pmin(lam, lambda[1:length(lam)])
                MAX <- min(sum(lambda), MAX)
        }

        Lp <- length(lambda)
        while(lambda[Lp] == 0) {
                Lp <- Lp-1
        }

        lambda <- lambda[1:Lp]

        f <- 1:(MAX+1)
        if(Lp - 1 > 1) {
                for(i in 2:(Lp-1)) {
                        for(j in ((i+1):(MAX+1))) { # aj pred tento for cyklus jeden "if" navyse?
                                f[j] <- f[j] + f[j-i]
                        }
                }
        }

        f <- f[length(f)]

        D  <- matrix(0, f, 1)
        Sx <- matrix(0,f, n)
        Sx[1,] <- 1

        xn <- matrix(1, n, MAX + 1)
        if(is.matrix(x) && ncol(x) > 1) { # tuto podmienku este skontrolovat!
                x <- t(x)
        }
        # if size(x,2) > 1
        # x = x.';
        # end

        for(i in 2:(MAX+1)) { # tento cyklus este skontrolovat!
                xn[,i] <- xn[,i-1] * x
        }
        # for i = 2:MAX+1
        # xn(:,i) = xn(:,i-1) .* x;
        # end
        prodx <- numeric()
        if(is.vector(x)) {
                prodx <- cumprod(x)
        } else if(is.matrix(x)) {
                prodx <- cumprod(x)
                dim(prodx) <- dim(x)
        }

        XY <- FALSE
        if(!missing(a) && !missing(b) && !missing(x) && !missing(y) && !missing(alpha) && (!missing(MAX) || !missing(lam)) && length(y) > 0) {
                XY <- TRUE
        }

        Sy <- numeric()
        yn <- numeric()
        prody <- numeric()
        if(XY) {
                Sy <- Sx
                yn <- matrix(1, n, MAX + 1)
                if(is.matrix(y) && ncol(y) > 1) { # tuto podmienku este skontrolovat!
                        y <- t(y)
                }
                for(i in 2:(MAX+1)) {
                        yn[,i] <- yn[,i-1] * y
                }
                prody <- numeric()
                if(is.vector(y)) {
                        prody <- cumprod(y)
                } else if(is.matrix(y)) {
                        prody <- cumprod(y)
                        dim(prody) <- dim(y)
                }
        }

        l <- rep(0, Lp)
        z <- matrix(1, na, Lp)
        kt <- -(1:Lp)
        cc1 <- 0
        ss <- matrix(0, na, MAX + 1)
        ss[,1] <- 1
        sl <- 1
        h <- 1
        ww <- rep(1, Lp)
        heap <- lambda[1] + 2
        d <- rep(0, Lp)

        while(h > 0) {

                if((l[h] < lambda[h]) && (h == 1 || l[h] < l[h-1]) && (MAX >= sl) && (z[h] != 0)) {

                        l[h] <- l[h] + 1

                        if(l[h]==1 && h>1 && h<n) {

                                D[ww[h]] <- heap
                                ww[h] <- heap
                                m <- min(lambda[h], MAX - sl + l[h])
                                heap <- heap + min(m, l[h-1])

                        } else {

                                ww[h] <- ww[h] + 1
                        }

                        w <- ww[h]
                        c <- -(h-1) / alpha + l[h] - 1

                        zn <- numeric()
                        dn <- numeric()
                        if(is.vector(a)) {
                                zn <- (a + c) * alpha
                        } else if(is.matrix(a)) {
                                zn <- apply(a + c, 1, prod) * alpha
                        }
                        if(is.vector(b)) {
                                dn <- (b + c) * (kt[h] + h + 1)
                        } else if(is.matrix(b)) {
                                dn <- apply(b + c, 1, prod) * (kt[h] + h + 1)
                        }

                        delta <- numeric()

                        if(XY) {
                                zn <- zn * alpha * l[h]
                                dn <- dn * (n + alpha * c)
                                if(h - 1 > 0) {
                                        for(j in 1:(h-1)) {
                                                delta <- kt[j] - kt[h]
                                                zn    = zn * delta
                                                dn    = dn * (delta - 1)
                                        }
                                }
                        }

                        kt[h] <- kt[h] + alpha

                        if(h - 1 > 0) {
                                for(j in 1:(h-1)) {
                                        delta <- kt[j] - kt[h]
                                        zn <- zn * delta
                                        dn <- dn * (delta + 1)
                                }
                        }

                        z[,h] <- z[,h] * zn / dn
                        sl <- sl + 1

                        if(h < n) {

                                if(h > 1) {
                                        d[h-1] <- d[h-1] - 1
                                }

                                d[h] <- l[h]
                                cc <- prod(h + 1 - alpha + kt[1:h]) / prod(h + kt[1:h])
                                pp <- l[1]
                                k <- 2

                                while(k <= h && l[k] > 1) {
                                        pp <- D[pp] + l[k] - 2
                                        k <- k + 1
                                }

                                Sx[w,h] <- cc * prodx[h] * Sx[pp,h]

                                if(XY) {
                                        Sy[w,h] <- cc * prody[h] * Sy[pp,h]
                                }

                                g <- which(d[1:h] > 0)
                                lg <- length(g)
                                slm <- 1
                                nhstrip <- prod(d[g]+1) - 1
                                mu <- l
                                mt <- kt
                                blm <- rep(1, lg)
                                lmd <- l[g] - d[g]

                                for(i in 1:nhstrip) {

                                        j <- lg
                                        gz <- g[lg]

                                        while(mu[gz] == lmd[j]) {

                                                mu[gz] <- l[gz]
                                                mt[gz] <- kt[gz]
                                                slm <- slm - d[gz]
                                                j <- j - 1
                                                gz <- g[j]
                                        }

                                        t <- kt[gz] - mt[gz]
                                        blm[j] <- blm[j] * (1 + t)
                                        dn <- t + alpha

                                        if(gz - 1 > 0) {
                                                for(r in 1:(gz-1)) {

                                                        q1 <- mt[r] - mt[gz]
                                                        q2 <- kt[r] - mt[gz]
                                                        blm[j] <- blm[j] * (q1 + alpha - 1) * (1 + q2)
                                                        dn <- dn * q1 * (alpha + q2)
                                                }
                                        }

                                        blm[j] <- blm[j] / dn
                                        mu[gz] <- mu[gz] - 1
                                        mt[gz] <- mt[gz] - alpha
                                        slm <- slm + 1

                                        if(j < lg) {

                                                blm[(j+1):length(blm)] <- blm[j]
                                        }

                                        nmu <- mu[1] + 1

                                        if(h-(mu[h]==0) > 1) {
                                                for(k in 2:(h-(mu[h]==0))) { # overit podmienku for cyklu

                                                        nmu <- D[nmu] + mu[k] - 1
                                                }
                                        }

                                        for(k in (h+1):n) {
                                                Sx[w,k] <- Sx[w,k] + blm[j] * Sx[nmu,k-1] * xn[k,slm]
                                        }

                                        if(XY) {

                                                for(k in (h+1):n) {
                                                        Sy[w,k] <- Sy[w,k] + blm[j] * Sy[nmu,k-1] * yn[k,slm]
                                                }
                                        }
                                }

                                for(k in (h+1):n) {
                                        Sx[w,k] <- Sx[w,k] + Sx[w,k-1]
                                }

                                if(XY) {

                                        for(k in (h+1):n) {
                                                Sy[w,k] <- Sy[w,k] + Sy[w,k-1]
                                        }

                                        ss[,sl] <- ss[,sl] + z[,h] * Sx[w,n] * Sy[w,n]
                                } else {

                                        ss[,sl] <- ss[,sl] + z[,h] * Sx[w,n]
                                }
                        } else {

                                pp <- l[1] + 1 - l[n]
                                k <- 2

                                while(k <= n-1 && l[k] > l[n]) {
                                        pp <- D[pp] + l[k] - 1 - l[n]
                                        k <- k + 1
                                }

                                k <- as.numeric((l[n]==1))
                                cc1 <- k + (1 - k) * cc1

                                if(XY) {

                                        cc1 <- cc1 * (prod(1 + kt[1:(n-1)] - kt[n]) *
                                                              (1 / alpha + l[n] - 1) /
                                                              (prod(alpha + kt[1:(n-1)] -
                                                                kt[n]) * l[n]))^2 * prodx[n] * prody[n]

                                        ss[,sl] <- ss[,sl] + z[,n] * cc1 * Sx[pp,n] * Sy[pp,n]

                                } else {

                                        cc1 <- cc1 * prod(1 + kt[1:(n-1)] - kt[n]) *
                                                prodx[n] * (1 / alpha + l[n] - 1) /
                                                (prod(alpha + kt[1:(n-1)] - kt[n]) * l[n])

                                        ss[,sl] <- ss[,sl] + z[,n] * cc1 * Sx[pp,n]
                                }
                        }

                        if(h < Lp) {

                                z[,h+1] <- z[,h]
                                ww[h+1] <- w
                                h <- h + 1
                        }

                } else {

                        sl <- sl - l[h]
                        l[h] <- 0
                        kt[h] <- -h
                        h <- h - 1
                }

        }

        s <- rowSums(ss)
        result <- list("s" = s, "ss" = ss)

        return(result)
}
