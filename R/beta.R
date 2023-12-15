######
## Beta distribution helpers
######


#' Plot the beta distribution with 95% area interval.
#'
#' @param a the first shape parameter
#' @param b the second shape parameter
#' @param level The leve of the interval.
#' @param ... Additional arguments to pass to \code{\link[graphics]{plot}}.
#'
#' @author David Gerard
#'
#' @return Plots the beta density with \code{level} interval shaded, along
#'     with the vertical line of the mean.
#'
#' @examples
#' plot_beta(20, 40)
#' plot_beta(1/5, 2/5, ylim = c(0, 5))
#'
#' @export
plot_beta <- function(a, b, level = 0.8, ...) {
  x <- seq(0, 1, length.out = 500)
  y <- stats::dbeta(x = x, shape1 = a, shape2 = b)

  l <- (1 - level) / 2
  ran <- stats::qbeta(p = c(l, 1 - l), shape1 = a, shape2 = b)

  xsub <- c(ran[[1]], x[x >= ran[[1]] & x <= ran[[2]]], ran[[2]], ran[[2]])
  ysub <- c(0, y[x >= ran[[1]] & x <= ran[[2]]], 0, 0)

  oldpar <- graphics::par(pch = 16,
                          mar = c(3, 3, 0.5, 0.5),
                          mgp = c(1.8, 0.4, 0),
                          tcl = -0.25)
  on.exit(graphics::par(oldpar), add = TRUE)

  graphics::plot(x = x, y = y, type = "l", xlab = expression(gamma), ylab = expression(plain(f)(gamma)), ...)
  graphics::abline(v = a / (a + b), lty = 2, col = 2)
  graphics::polygon(x = xsub, y = ysub, col = "#0000FF22", border = NA)
  graphics::grid()

  invisible(NULL)
}
