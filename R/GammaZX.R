#' @title
#' Gamma function valid in the entire complex plane
#'
#' @description
#' \code{GammaZX(z, funmode)} evaluates Gamma function valid in the entire complex plane, the argument \code{z}
#' may be complex and of any size.
#'
#' @family Utility Function
#'
#' @param z complex argument.
#' @param funmode function mode, \code{funmode = 0} for \eqn{ln[gamma(z)]} and \code{funmode = 1} for \eqn{gamma(z)}.
#'
#' @return  Function returns \eqn{ln[gamma(z)]} or \eqn{gamma(z)}.
#'
#' @note Ver.: 01-Oct-2018 13:00:11 (consistent with Matlab CharFunTool v1.3.0, 20-Jul-2018 16:13:31).
#'
#' @export
#'
GammaZX <- function(z, funmode) {
        ## CHECK THE INPUT PARAMETERS
        if(missing(funmode)) {
                funmode <- 1
        }

        ## ALGORITHM
        sz <- dim(z)
        z <- c(z)

        x <- Re(z)
        y <- Im(z)
        x1 <- 0
        pi <- 3.141592653589793
        a <- c(8.333333333333333e-2, -2.777777777777778e-3, 7.936507936507937e-4, -5.952380952380952e-4,
               8.417508417508418e-4, -1.917526917526918e-3, 6.410256410256410e-3, -2.955065359477124e-2,
               1.796443723688307e-1, -1.39243221690590)

        if ((all(y == 0 && x == as.integer(x))) && all(x <= 0)) {
                gr <- 1e300
                gi <- 0
                g  <- gr + 1i*gi
                return(g)
        } else if (all(x < 0)) {
                x1 <- x
                x <- -x
                y <- -y
        }

        x0 <- x
        if(all(x <= 7)) {
                na <- as.integer(7 - x)
                x0 <- x + na
        }

        z1 <- sqrt(x0 * x0 + y * y)
        th <- atan(y / x0)
        gr <- (x0 - 0.5) * log(z1) - th * y - x0 + 0.5 * log(2 * pi)
        gi <- th * (x0 - 0.5) + y * log(z1) - y
        for(k in 1:10) {
                t <- z1 ^ (1 - 2 * k)
                gr <- gr + a[k] * t * cos((2 * k - 1) * th)
                gi <- gi - a[k] * t * sin((2 * k - 1) * th)
        }

        if(all(x <= 7)) {
                gr1 <- 0
                gi1 <- 0
                for(j in 0:(na-1)) {
                        gr1 <- gr1 + 0.5 * log((x + j) ^ 2 + y * y)
                        gi1 <- gi1 + atan(y / (x + j))
                }

                gr <- gr - gr1
                gi <- gi - gi1
        }

        if(all(x1 < 0)) {
                z1 <- sqrt(x * x + y * y)
                th1 <- atan(y / x)
                sr <- -sin(pi * x) * cosh(pi * y)
                si <- -cos(pi * x) * sinh(pi * y)
                z2 <- sqrt(sr * sr + si * si)
                th2 <- atan(si / sr)
                if(any(sr < 0)) {
                        th2 <- pi + th2
                }
                gr = log(pi / (z1 * z2)) - gr
                gi = -th1 - th2 - gi
        }

        if(funmode == 1) {
                g0 <- exp(gr)
                gr <- g0 * cos(gi)
                gi <- g0 * sin(gi)
        }
        g <- gr + 1i * gi

        if(all(Im(z) == 0)) {
                g <- Re(g)
        }

        dim(g) <- sz

        return(g)
}

