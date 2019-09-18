#' @title
#' Gauss hypergeometric function 2F1(a,b,c,z)
#'
#' @description
#' \code{Hypergeom2F1(a, b, c, z)} computes the Gauss hypergeometric function \eqn{2F1(a,b,c,z)},
#' for the real parameters \code{a}, \code{b} and \code{c} (here assumed to be scalars),
#' and the complex argument \code{z} (could be scalar, vector or array).
#'
#' @details
#' The functions based on \code{Hypergeom2F1} is badly conditioned for when \eqn{c} is negative integer.
#' In such situations, approximate the function value by using the noninteger parameter \eqn{c},
#' say \eqn{c = c + eps}, for some small \eqn{eps}.
#'
#' @references
#' The algorithm is based on a Fortran program in
#'
#' S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
#' Converted by Ben Barrowes (barrowes@alum.mit.edu).
#'
#' @param a parameter - real scalar.
#' @param b parameter - real scalar.
#' @param c parameter - real scalar.
#' @param z complex argument - scalar, vector or array.
#'
#' @family Utility Function
#'
#' @return
#' The Gauss hypergeometric function \eqn{2F1(a,b,c,z)}.
#'
#' @note Ver.: 06-Oct-2018 18:36:05 (consistent with Matlab CharFunTool v1.3.0, 18-Aug-2018 18:32:27).
#'
#' @example R/Examples/example_Hypergeom2F1.R
#'
#' @export
#'
Hypergeom2F1 <- function(a, b, c, z) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(a) || missing(b) || missing(c) || missing(z)) {
                stop("All input paramaters must be entered.")
        }

        l_max <- max(c(length(a), length(b), length(c)))
        if (l_max > 1) {
                if (length(b) == 1) {
                        b <- rep(b, l_max)
                }
                if (length(a) == 1) {
                        a <- rep(a, l_max)
                }
                if (length(c) == 1) {
                        c <- rep(c, l_max)
                }
                if ((any(lengths(list(a, b, c)) < l_max))) {
                        stop("Input size mismatch.")
                }
        }

        na <- length(a)

        sz <- dim(z)
        z <- c(z)
        nz <- length(z)

        if(nz >= 1 && na == 1) {
                f <- rep(0, nz)
                for(i in 1:nz) {
                        f[i] <- hygfz(a, b, c, z[i])
                }
                dim(f) <- sz
                return(f)
        } else if(nz == 1 && na > 1) {
                f <- rep(0, na)
                for(i in 1:na) {
                        f[i] <- hygfz(a[i], b[i], c[i], z)
                }
                return(f)
        } else {
                stop("Input size mismatch.")
        }
}

hygfz <- function(a, b, c, z) {
        #  HYGFZ Compute the hypergeometric function for a complex argument,
        #  F(a,b,c,z)
        #  Input:
                # a --- Parameter
                # b --- Parameter
                # c --- Parameter,  c <> 0,-1,-2,...
                # z --- Complex argument
        # Output:
                # ZHF --- F(a,b,c,z)

        ## ALGORITHM
        zw <- 0
        x <- Re(z)
        y <- Im(z)
        eps <- 1e-15
        l0 <- (c == round_tw0(c) && Re(c) < 0)
        l1 <- (abs(1 - x) < eps && y == 0 && Re(c - a - b) <= 0)
        l2 <- (abs(z + 1) < eps && abs(c - a + b - 1) < eps)
        l3 <- (a == round_tw0(a) && Re(a) < 0)
        l4 <- (b == round_tw0(b) && Re(b) < 0)
        l5 <- (c - a == round_tw0(c - a) && Re(c - a) <= 0)
        l6 <- (c - b == round_tw0(c - b) && Re(c - b) <= 0)
        aa <- a
        bb <- b
        a0 <- abs(z)
        if(a0 > 0.95) {
                eps <- 1e-8
        }

        pi <- 3.141592653589793
        el <- 0.5772156649015329
        if(l0 || l1) {

                cat(1,'%s \n','the hypergeometric series is divergent')
                return()

        }

        if(a0 == 0 || a == 0 || b == 0) {

                zhf <- 1

        } else if(z == 1 && Re(c - a - b) > 0) {

                gc <- gamma(c)
                gcab <- gamma(c - a - b)
                gca <- gamma(c - a)
                gcb <- gamma(c - b)
                zhf <- gc * gcab / (gca * gcb)

        } else if(l2) {

                g0 <- sqrt(pi) * 2^(-a)
                g1 <- gamma(c)
                g2 <- gamma(1 + a/2 - b)
                g3 <- gamma(0.5 + 0.5 * a)
                zhf <- g0 * g1 / (g2 * g3)

        } else if(l3 || l4) {

                if(l3) {
                        nm <- round_tw0(abs(a))
                }

                if(l4) {
                        nm <- round_tw0(abs(b))
                }

                zhf <- 1
                zr <- 1
                for(k in 1:nm) {
                        zr <- zr * (a + k - 1) * (b + k - 1) / (k * (c + k - 1)) * z
                        zhf <- zhf + zr
                }

        } else if(l5 || l6) {

                if(l5) {
                        nm <- round_tw0(abs(c - a))
                }

                if(l6){
                        nm <- round_tw0(abs(c - b))
                }

                zhf <- 1
                zr <- 1
                for(k in 1:nm) {
                        zr <- zr * (c - a + k - 1) * (c - b + k - 1) / (k * (c + k - 1)) * z
                        zhf <- zhf + zr
                }

                zhf  = (1 - z) ^ (c - a - b) * zhf

        } else if(a0 <= 1) {

                if (x < 0) {

                        z1 <- z / (z - 1)

                        if (c > a && b < a && b > 0) {
                                a <- bb
                                b <- aa
                        }

                        zc0 <- 1 / ((1 - z) ^ a)
                        zhf <- 1
                        zr0 <- 1
                        for(k in 1:500) {
                                zr0 <- zr0 * (a + k - 1) * (c - b + k - 1) / (k * (c + k - 1)) * z1
                                zhf <- zhf + zr0
                                if (abs(zhf - zw) < abs(zhf) *eps) {
                                        break
                                }
                                zw <- zhf
                        }
                        zhf <- zc0 * zhf

                } else if(a0 >= 0.9) {

                        gm <- 0
                        mcab <- round_tw0(c - a - b + eps * sign(c - a - b))

                        if (abs(c - a - b - mcab) < eps) {

                                m <- round_tw0(c - a - b)
                                ga <- gamma(a)
                                gb <- gamma(b)
                                gc <- gamma(c)
                                gam <- gamma(a + m)
                                gbm <- gamma(b + m)
                                pa <- psigamma(a)
                                pb <- psigamma(b)
                                if(m != 0) {
                                        gm <- 1
                                }

                                for(j in 1:(abs(m)-1)) {
                                        gm <- gm * j
                                }

                                rm <- 1

                                for(j in 1:(abs(m))) {
                                        rm <- rm * j
                                }

                                zf0 <- 1
                                zr0 <- 1
                                zr1 <- 1
                                sp0 <- 0
                                sp <- 0

                                if (m >= 0) {

                                        zc0 <- gm * gc / (gam * gbm)
                                        zc1 <- -gc * (z - 1) ^ m / (ga * gb * rm)
                                        for(k in 1:(m-1)) {
                                                zr0 <- zr0 * (a + k - 1) * (b + k - 1) / (k * (k - m)) * (1 - z)
                                                zf0 <- zf0 + zr0
                                        }

                                        for(k in 1:m) {
                                                sp0 <- sp0 + 1 / (a + k - 1) + 1 / (b + k - 1) - 1 / k
                                        }

                                        zf1 <- pa + pb + sp0 + 2 * el + log(1 - z)

                                        for(k in 1:500) {
                                                sp <- sp + (1 - a) / (k * (a + k - 1)) + (1 - b) / (k * (b + k - 1))
                                                sm <- 0

                                                for(j in 1:m) {
                                                        sm <- sm + (1 - a) / ((j + k) * (a + j + k - 1)) + 1 / (b + j + k - 1)
                                                }

                                                zp <- pa + pb + 2 * el + sp + sm + log(1 - z)
                                                zr1 <- zr1 * (a + m + k - 1) * (b + m + k - 1) / (k * (m + k)) * (1 - z)
                                                zf1 <- zf1 + zr1 * zp

                                                if(abs(zf1 - zw) < abs(zf1) * eps) {
                                                        break
                                                }

                                                zw <- zf1
                                        }

                                        zhf <- zf0 * zc0 + zf1 * zc1

                                } else if(m < 0) {

                                        m  <- -m
                                        zc0 <- gm * gc / (ga * gb * (1 - z) ^ m)
                                        zc1 <-  -(-1) ^ m * gc / (gam * gbm * rm)

                                        for(k in 1:(m-1)) {
                                                zr0 <- zr0 * (a - m + k - 1) * (b - m + k - 1) / (k * (k-m)) * (1 - z)
                                                zf0 <- zf0 + zr0
                                        }

                                        for(k in 1:m) {
                                                sp0 <- sp0 + 1 / k
                                        }

                                        zf1 <- pa + pb - sp0 + 2 * el + log(1 - z)
                                        for(k in 1:500) {
                                                sp <- sp + (1 - a) / (k * (a + k - 1)) + (1 - b) / (k * (b + k - 1))
                                                sm <- 0

                                                for(j in 1:m) {
                                                        sm <- sm + 1 / (j + k)
                                                }

                                                zp <- pa + pb + 2 * el + sp - sm + log(1 - z)
                                                zr1 <- zr1 * (a + k - 1) * (b + k - 1) / (k * (m + k)) * (1 - z)
                                                zf1 <- zf1 + zr1 * zp

                                                if(abs(zf1 - zw) < abs(zf1) * eps) {
                                                        break
                                                }

                                                zw <- zf1
                                        }

                                        zhf <- zf0 * zc0 + zf1 * zc1
                                }

                        } else {

                                ga <- gamma(a)
                                gb <- gamma(b)
                                gc <- gamma(c)
                                gca <- gamma(c - a)
                                gcb <- gamma(c - b)
                                gcab <- gamma(c - a - b)
                                gabc <- gamma(a + b - c)
                                zc0 <- gc * gcab / (gca * gcb)
                                zc1 <- gc * gabc / (ga * gb) * (1-z)^(c-a-b)
                                zhf <- 0
                                zr0 <- zc0
                                zr1 <- zc1
                                for(k in 1:500) {
                                        zr0 <- zr0 * (a + k - 1) * (b + k - 1) / (k * (a + b - c + k)) * (1 - z)
                                        zr1 <- zr1 * (c - a + k - 1) * (c - b + k - 1) / (k * (c - a - b + k)) * (1 - z)
                                        zhf <- zhf + zr0 + zr1

                                        if (abs(zhf - zw) < abs(zhf) * eps) {
                                                break
                                        }

                                        zw <- zhf
                                }

                                zhf <- zhf + zc0 + zc1
                        }

                } else {

                        z00 <- 1

                        if (Re(c - a) < Re(a) && Re(c - b) < Re(b)) {

                                z00 <- (1 - z) ^ (c - a - b)
                                a <- c - a
                                b <- c - b
                        }

                        zhf <- 1
                        zr <- 1

                        for(k in 1:1500) {
                                zr <- zr * (a + k - 1) * (b + k - 1) / (k * (c + k - 1)) * z
                                zhf <- zhf + zr

                                if(abs(zhf - zw) <= abs(zhf) * eps) {
                                        break
                                }

                                zw <- zhf
                        }

                        zhf <- z00 * zhf

                }

        } else if(a0 > 1) {

                mab <- round_tw0(a - b + eps * sign(a - b))
                if(abs(a - b - mab) < eps && a0 <= 1.1) {

                        b <- b + eps

                }

                if (abs(a - b - mab) > eps) {

                        ga <- gamma(a)
                        gb <- gamma(b)
                        gc <- gamma(c)
                        gab <- gamma(a - b)
                        gba <- gamma(b - a)
                        gca <- gamma(c - a)
                        gcb <- gamma(c - b)
                        zc0 <- gc * gba / (gca * gb * (-z) ^a)
                        zc1 <- gc * gab / (gcb * ga *(-z) ^b)
                        zr0 <- zc0
                        zr1 <- zc1
                        zhf <- 0

                        for(k in 1:500) {
                                zr0 <- zr0 * (a + k - 1) * (a - c + k) / ((a - b + k) * k *z)
                                zr1 <- zr1 * (b + k - 1) * (b - c + k) / ((b - a + k) * k *z)
                                zhf <- zhf + zr0 + zr1
                                if (abs((zhf - zw) / zhf) <= eps) {
                                        break
                                }

                                zw <- zhf
                        }

                        zhf <- zhf + zc0 + zc1

                } else {

                        if (a - b < 0) {

                                a <- bb
                                b <- aa

                        }

                        ca <- c - a
                        cb <- c - b
                        nca <- round_tw0(ca + eps * sign(ca))
                        ncb <- round_tw0(cb + eps * sign(cb))

                        if(abs(ca - nca) < eps || abs(cb - ncb) < eps) {
                                c <- c + eps
                        }

                        ga <- gamma(a)
                        gc <- gamma(c)
                        gcb <- gamma(c - b)
                        pa <- digamma(a)
                        pca <- digamma(c - a)
                        pac <- digamma(a - c)
                        mab <- round_tw0(a - b + eps)
                        zc0 <- gc / (ga * (-z) ^ b)
                        gm <- gamma(a - b)
                        zf0 <- gm / gcb * zc0
                        zr <- zc0

                        for(k in 1:(mab-1)) {
                                zr <- zr * (b + k - 1) / (k * z)
                                t0 <- a - b - k
                                g0 <- gamma(t0)
                                gcbk <- gamma(c - b - k)
                                zf0 <- zf0 + zr * g0 / gcbk
                        }

                        if(mab == 0) {
                                zf0 <- 0
                        }

                        zc1 <- gc / (ga * gcb * (-z) ^a)
                        sp <- -2 * el - pa - pca
                        for(j in 1:mab) {
                                sp <- sp + 1 / j
                        }

                        zp0 <- sp + log(-z)
                        sq <- 1
                        for(j in 1:mab) {
                                sq <- sq * (b + j - 1) * (b - c + j) / j
                        }

                        zf1 <- (sq * zp0) * zc1
                        zr <- zc1
                        rk1 <- 1
                        sj1 <- 0
                        # VW added next line: w0 = 0
                        w0 <- 0

                        for(k in 1:10000) {
                                zr <- zr / z
                                rk1 <- rk1 * (b + k - 1) * (b - c + k) / (k * k)
                                rk2 <- rk1

                                for(j in (k+1):(k+mab)) {
                                        rk2 <- rk2 * (b + j - 1) * (b - c + j) / j
                                }

                                sj1 <- sj1 + (a - 1) / (k * (a + k - 1)) + (a - c - 1) / (k * (a - c + k - 1))
                                sj2 <- sj1

                                for(j in (k+1):(k+mab)) {
                                        sj2 <- sj2 + 1 / j
                                }

                                zp <- -2 * el - pa - pac + sj2 - 1 / (k + a - c) - pi / tan(pi * (k + a - c)) + log(-z)
                                zf1 <- zf1 + rk2 * zr * zp
                                ws <- abs(zf1)

                                if (abs((ws - w0) / ws) < eps) {
                                        break
                                }

                                w0 <- ws
                        }

                        zhf <- zf0 + zf1
                }

        }

        return(zhf)
}

# auxiliary function, an equivalent to Matlab fix function (round also complex numbers towards zero)
round_tw0 <- function(z) {

        result <- as.integer(Re(z))+as.integer(Im(z))*1i

        return(result)
}
