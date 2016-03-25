#' Constructs a block diagonal matrix
#'
#' The function takes a \code{list} and constructs a block diagonal matrix with
#' the elements of the list on the diagonal. If \code{d} is not a list then
#' \code{d} will be repeated \code{n} times and written on the diagonal (a wrapper for \code{kronecker()})
#'
#' @param d a \code{list} of matrices or vectors, or a matrix or vector
#' @param n number of repetitions
#'
#' @return returns a matrix with the elements of the list or the repetitions of the supplied matrix or vector on the diagonal.
#'
#' @author Wolfgang Lederer, \email{wolfgang.lederer@gmail.com}
#'
#' @seealso \code{\link[base]{diag}}, \code{\link[base]{kronecker}}
#'
#' @examples
#' a <- matrix(rep(1, 4), nrow = 2)
#' b <- matrix(rep(2, 6), nrow = 2)
#' e <- c(3, 3, 3, 3)
#' f <- t(e)
#' d <- list(a, b, e, f)
#' diag.block(d)
#' diag.block(a, 3)
#'
#' @export




diag.block <- function(d, n) {
  if (is.list(d)) {
    d.row <- sapply(d, NROW)
    d.col <- sapply(d, NCOL)
    d.diag <- matrix(0, nrow = sum(d.row), ncol = sum(d.col))
    d.row <- c(0, cumsum(d.row))
    d.col <- c(0, cumsum(d.col))
    for (i in 1:length(d)) {
      d.diag[(d.row[i] + 1):d.row[i + 1], (d.col[i] + 1):d.col[i + 1]] <- as.matrix(d[[i]])
    }
  }
  if (!is.list(d)) {
    d.diag <- kronecker(diag(1, n), d)
  }
  return(d.diag)
}

