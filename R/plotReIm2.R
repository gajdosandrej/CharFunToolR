#' @title Plot real and imaginary part of f(x)
#'
#' @description
#' \code{plotReIm2(f, x)} plots values of complex functions \code{f} in real arguments \code{x} or in ranges \code{[xmin, xmax]} with steps \code{dx}.
#'
#' @importFrom graphics plot lines legend grid
#' @importFrom randomcoloR distinctColorPalette
#'
#' @param f list of functions.
#' @param x list of numerical values (numerical vectors).
#' @param title an overall title for the plot.
#' @param xmin vector of numbers, default values \code{xmin = -1}.
#' @param xmax vector of numbers, default values \code{xmax = 1}.
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
#' @note Ver.: 22-Sep-2018 20:41:29.
#'
#' @export
#'
plotReIm2 <- function(f, x, title = "", xmin, xmax, dx = 0.1, xlab = "t", ylab = "f(t)",
                      px = 0, py = 0, labelx = TRUE, labely = TRUE) {
        if(missing(xmin)) {
                xmin <- vector()
        }
        if(missing(xmax)) {
                xmax <- vector()
        }

        x_new <- list()
        f_real <- list()
        f_imag <- list()
        if(!missing(x)) {
            x_new <- x
        }
        xmin_new <- xmin
        xmax_new <- xmax

        for(i in 1:length(f)) {
                if(length(x_new) < i) {
                        if (length(xmin) < i && length(xmax) < i) {
                                xmin_new[i] <- -1
                                xmax_new[i] <- 1
                        } else if (length(xmin) < i) {
                                xmin_new[i] <- xmax - 2
                        } else if (length(xmin) < i) {
                                xmax_new[i] <- xmin + 2
                        }
                        x_new[[i]] <- seq(xmin_new[i], xmax_new[i], by = dx)
                } else {
                        if (length(xmin) < i) {
                                xmin_new[i] <- min(x[[i]])
                        }
                        if (length(xmax) < i) {
                                xmax_new[i] <- max(x[[i]])
                        }
                        x_new[[i]] <- x[[i]][x[[i]] >= xmin_new[i] && x[[i]] <= xmax_new[i]]
                }
                f_new <- f[[i]](x_new[[i]])
                f_real[[i]] <- Re(f_new)
                f_imag[[i]] <- Im(f_new)
        }

                plot(NULL, main = title, xlab = xlab, ylab = ylab, xlim = c(min(xmin_new), max(xmax_new)),
                     ylim = c(min(c(unlist(f_real), unlist(f_imag))), max(c(unlist(f_real), unlist(f_imag)))))

                if (labelx) {
                        for (xi in c(px)) {
                                lines(c(xi, xi),
                                        c(min(c(unlist(f_real), unlist(f_imag))) - 1, max(c(unlist(f_real), unlist(f_imag))) + 1),
                                        lty = 1,
                                        lwd = 2,
                                        col = "gray"
                                )
                        }
                }
                if (labely) {
                        for (yi in c(py)) {
                                lines(
                                        c(min(xmin_new), max(xmax_new)),
                                        c(yi, yi),
                                        lty = 1,
                                        lwd = 2,
                                        col = "gray"
                                )
                        }
                }

                lines_col <- randomcoloR::distinctColorPalette(2 * length(f))

                for(j in 1:length(f)) {
                        lines(x_new[[j]], f_real[[j]], lty = 1, lwd = 2, col = lines_col[2 * j - 1])
                        lines(x_new[[j]], f_imag[[j]], lty = 1, lwd = 2, col = lines_col[2 * j])
                }

                legend("bottomright",
                       legend = c("Re(f))",rep("", length(f) - 1),"Im(f)",rep("", length(f) - 1)),
                       fill = lines_col,
                       bg = "white", ncol = 2, cex = 0.75)
                grid()

        }
