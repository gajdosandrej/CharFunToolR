#' @title Plot real and imaginary part of f(x)
#'
#' @description
#' \code{plotReIm(f, x, xmin, xmax)} plots values of complex function \code{f} in real arguments \code{x} or in range \code{[xmin, xmax]} with step \code{dx}.
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
#' @param px numerical values (number or vector) in which are auxiliary lines on axis \eqn{x} rendered, default value \code{px = 0}.
#' @param py numerical values (number or vector) in which are auxiliary lines on axis \eqn{y} rendered, default value \code{py = 0}.
#' @param labelx rendering auxiliary line, default value \code{TRUE}.
#' @param labely rendering auxiliary line, default value \code{TRUE}.
#'
#' @details
#' Function with optional parameters must be pass as lambda function \code{f = function(x) f(x, optional parameters)}.
#'
#'
#' @return Graf of function \code{f} in argument \code{x} or in range \code{[xmin, xmax]} with step \code{dx}.
#'
#' @family Graphic Function
#'
#' @note Ver.: 22-Sep-2018 20:32:20.
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
                 labely = TRUE
                 ) {
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
                        x <- x[x >= xmin & x <= xmax]
                }

                if (missing(title)) {
                        title <- "Function f"
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

                plot(NULL, main = title, xlab = xlab, ylab = ylab, xlim = c(xmin, xmax), ylim = c(min(real, imag), max(real, imag)))

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

                lines(x,
                      real,
                      lty = 1,
                      lwd = 2,
                      col = "blue")

                lines(x,
                      imag,
                      lty = 1,
                      lwd = 2,
                      col = "red")

                legend("bottomright",
                       c("Re(f))", "Im(f)"),
                       fill = c("blue", "red"),
                       cex = 0.75)
                grid()

        }
