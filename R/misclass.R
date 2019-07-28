#' Generates misclassified data
#'
#' Takes a \code{data.frame} and produces misclassified data.
#' Probabilities for the missclassification are given in \code{mc.matrix}.
#'
#' @param data.org \code{data.frame} containing the factor variables. Must be \code{factor}s.
#' @param mc.matrix a \code{list} of matrices giving the probabilities for the
#' misclassification. Names of the \code{list} must correspond to the variable
#' names in \code{data.org}. The \code{colnames} must be named according to the factor levels.
#' @param k the exponent for the misclassification matrix
#'
#' @return A \code{data.frame} containing the misclassified variables
#'
#' @author Wolfgang Lederer, \email{wolfgang.lederer@gmail.com}
#'
#' @seealso \code{\link[simex]{mcsimex}}, \code{\link[simex]{mc.matrix}}, \code{\link[simex]{diag.block}}
#'
#' @examples
#'
#' x1 <- factor(rbinom(100, 1, 0.5))
#' x2 <- factor(rbinom(100, 2, 0.5))
#'
#' p1 <- matrix(c(1, 0, 0, 1), nrow = 2)
#' p2 <- matrix(c(0.8, 0.1, 0.1, 0.1, 0.8, 0.1, 0.1, 0.1, 0.8), nrow = 3)
#'
#' colnames(p1) <- levels(x1)
#' colnames(p2) <- levels(x2)
#'
#' x <- data.frame(x1 = x1, x2 = x2)
#' mc.matrix <- list(x1 = p1, x2 = p2)
#'
#' x.mc <- misclass(data.org = x, mc.matrix = mc.matrix, k = 1)
#'
#' identical(x[, 1], x.mc[, 1]) # TRUE
#' identical(x[, 2], x.mc[, 2]) # FALSE
#'
#'
#' @export
#'

misclass <- function(data.org, mc.matrix, k = 1) {
  if (!is.list(mc.matrix))
    stop("mc.matrix must be a list", call. = FALSE)
  if (!is.data.frame(data.org))
    stop("data.org must be a dataframe", call. = FALSE)
  if (!all(names(mc.matrix) %in% colnames(data.org)))
    stop("Names of mc.matrix and colnames of data.org do not match",
         call. = FALSE)
  if (k < 0)
    stop("k must be positive")
  data.mc <- data.org
  factors <- lapply(data.org, levels)
  ev <- lapply(mc.matrix, eigen)
  data.names <- colnames(data.org)
  for (j in data.names) {
    evalue <- ev[[c(j, "values")]]
    evectors <- ev[[c(j, "vectors")]]
    d <- diag(evalue)
    mc <- zapsmall(evectors %*% d^k %*% solve(evectors))
    dimnames(mc) <- dimnames(mc.matrix[[j]])
    for (i in factors[[j]]) {
      data.mc[[j]][data.org[[j]] == i] <- sample(x = factors[[j]],
                                                 size = length(data.org[[j]][data.org[[j]] == i]),
                                                 prob = mc[, i],
                                                 replace = TRUE)
    }
  }
  return(data.mc)
}
