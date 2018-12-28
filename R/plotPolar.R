#' @title Plot a line in polar coordinates
#'
#' @description
#'  \code{plotPolar(angle, radius, title, type, color)} plots a line in polar coordinates,
#'  with \code{angle} indicating the angle in radians and \code{radius} indicating the radius value for each point.
#'
#' @importFrom graphics plot lines par
#'
#' @param angle angle.
#' @param radius amplitude.
#' @param title an overall title for the plot.
#' @param type type of lines, dots.
#' @param color color.
#'
#' @family Graphic Function
#'
#' @note Ver.: 22-Sep-2018 20:35:05.
#'
#' @export
#'
plotPolar <-
  function(angle,
           radius,
           title = "",
           type = "l",
           color = "blue") {
    x <- radius * cos(angle)
    y <- radius * sin(angle)

    m <- max(abs(radius * exp(1i * angle)))
    m <- ceiling(m)

    par(pty = "s")

    plot(
      x = x,
      y = y,
      main = title,
      xlab = "",
      ylab = "",
      type = "l",
      col = "blue",
      xlim = c(-m, m),
      ylim = c(-m, m),
      asp = 1,
      lwd = 2
    )

    n <- 12

    for (i in seq(0, n / 2)) {
      a <- 2 * pi / n
      lines(
        m * c(cos(i * a),-cos(i * a)),
        m * c(sin(i * a),-sin(i * a)),
        lty = 1,
        lwd = 2,
        col = "gray"
      )
    }

    a <- seq(0, 2 * pi, length.out = 101)
    p <- 4
    m <- m / p

    for (i in seq(1:p)) {
      x <- m * i * cos(a)
      y <- m * i * sin(a)

      lines(
        x = x,
        y = y,
        type = "l",
        lwd = 2,
        col = "gray"
      )
    }

    par(pty = "m")
  }
