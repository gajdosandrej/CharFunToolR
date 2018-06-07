#' @title Plot real and imaginary part of f(x)
#'
#' @description
#' \code{plotReIm(f, x, xmin, xmax)} plots graf of complex function \code{f} in real arguments \code{x} or in range \code{[xmin, xmax]} with step \code{dx}.
#'
#' @importFrom graphics plot lines legend grid
#'
#' @param f function.
#' @param x numerical values (number or vector).
#' @param title an overall title for the plot.
#' @param xmin number, default value \code{xmin = -1}.
#' @param xmax number, default value \code{xmax = 1}.
#' @param dx positive number, default value \code{dx = 0.1}.
#' @param xlab a title for the \eqn{x} axis: see \link[graphics]{title}.
#' @param ylab a title for the \eqn{y} axis: see \link[graphics]{title}.
#' @param px numerical values (number or vector) in which are auxiliary line on axis \eqn{x} render default value \code{px = 0}.
#' @param py numerical values (number or vector) in which are auxiliary line on axis \eqn{y} render default value \code{py = 0}.
#' @param labelx rendering auxiliary line default value TRUE.
#' @param labely rendering auxiliary line default value TRUE.
#'
#' @details
#' Function with optional parameters must be pass as lambda function \code{f = function(x) f(x, optional parameters)}.
#'
#'
#' @return Graf of function \code{f} in argument \code{x} or in range \code{[xmin, xmax]} with step \code{dx}.
#'
#' @family Graphic Function
#'
#' @export
#'
plotReIm <-
        function(f,
                 x,
                 title,
                 xmin,
                 xmax,
                 dx = 0.1,
                 xlab,
                 ylab,
                 px = 0,
                 py = 0,
                 labelx = TRUE,
                 labely = TRUE) {
                if (missing(x)) {
                        if (missing(xmin) && missing(xmax)) {
                                xmin = -1
                                xmax = 1
                        } else if (missing(xmin)) {
                                xmin <- xmax - 2
                        } else if (missing(xmin)) {
                                xmax <- xmin + 2
                        }
                        x <- seq(xmin, xmax, by = dx)
                } else {
                        if (missing(xmin)) {
                                xmin <- min(x)
                        }
                        if (missing(xmax)) {
                                xmax <- max(x)
                        }
                        x <- x[x >= xmin && x <= xmax]
                }

                if (missing(title)) {
                        title <- "Function F(x)"
                }
                if (missing(xlab)) {
                        xlab <- "t"
                }
                if (missing(ylab)) {
                        ylab <- "f(t)"
                }

                f <- f(x)
                real <- Re(f)
                imag <- Im(f)

                plot(
                        x,
                        real,
                        main = title,
                        xlab = xlab,
                        ylab = ylab,
                        type = "l",
                        lwd = 2,
                        col = "blue",
                        xlim = c(xmin, xmax),
                        ylim = c(min(real, imag), max(real, imag))
                )
                lines(x,
                      imag,
                      lty = 1,
                      lwd = 2,
                      col = "red")

                if (labelx) {
                        for (xi in c(px)) {
                                lines(
                                        c(xi, xi),
                                        c(min(real, imag) - 1, max(real, imag) + 1),
                                        lty = 1,
                                        lwd = 2,
                                        col = "gray"
                                )
                        }
                }
                if (labely) {
                        for (yi in c(py)) {
                                lines(
                                        c(xmin, xmax),
                                        c(yi, yi),
                                        lty = 1,
                                        lwd = 2,
                                        col = "gray"
                                )
                        }
                }

                legend("bottomright",
                       c("Re(f))", "Im(f)"),
                       fill = c("blue", "red"))
                grid()

        }
