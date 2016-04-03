#' Misclassification in models using MCSIMEX
#'
#' Implementation of the misclassification MCSIMEX algorithm as described by Küchenhoff, Mwalili and Lesaffre.
#'
#' @aliases mcsimex print.mcsimex summary.mcsimex print.summary.mcsimex plot.mcsimex predict.mcsimex refit.mcsimex
#'
#' @param model the naive model, the misclassified variable must be a factor
#' @param SIMEXvariable vector of names of the variables for which the MCSIMEX-method should be applied
#' @param mc.matrix if one variable is misclassified it can be a matrix. If more
#' than one variable is misclassified it must be a list of the misclassification
#' matrices, names must match with the SIMEXvariable names, column- and row-names
#' must match with the factor levels. If a special misclassification is desired,
#' the name of a function can be specified (see details)
#' @param lambda vector of exponents for the misclassification matrix (without 0)
#' @param B number of iterations for each lambda
#' @param fitting.method \code{linear}, \code{quadratic} and \code{loglinear}
#' are implemented (first 4 letters are enough)
#' @param jackknife.estimation specifying the extrapolation method for jackknife
#' variance estimation. Can be set to \code{FALSE} if it should not be performed
#' @param asymptotic logical, indicating if asymptotic variance estimation should
#' be done, the option \code{x = TRUE} must be enabled in the naive model
#' @param x object of class 'mcsimex'
#' @param digits number of digits to be printed
#' @param object object of class 'mcsimex'
#' @param xlab optional name for the X-Axis
#' @param ylab vector containing the names for the Y-Axis
#' @param ask ogical. If \code{TRUE}, the user is asked for input, before a new figure is drawn
#' @param show vector of logicals indicating for which variables a plot should be produced
#' @param newdata optionally, a data frame in which to look for variables with
#' which to predict. If omitted, the fitted linear predictors are used
#' @param \dots arguments passed to other functions
#'
#' @return An object of class 'mcsimex' which contains:
#'
#' \item{coefficients}{corrected coefficients of the MCSIMEX model,}
#' \item{SIMEX.estimates}{the MCSIMEX-estimates of the coefficients for each lambda,}
#' \item{lambda}{the values of lambda,}
#' \item{model}{the naive model,}
#' \item{mc.matrix}{the misclassification matrix,}
#' \item{B}{the number of iterations,}
#' \item{extrapolation}{the model object of the extrapolation step,}
#' \item{fitting.method}{the fitting method used in the extrapolation step,}
#' \item{SIMEXvariable}{name of the SIMEXvariables,}
#' \item{call}{the function call,}
#' \item{variance.jackknife}{the jackknife variance estimates,}
#' \item{extrapolation.variance}{the model object of the variance extrapolation,}
#' \item{variance.jackknife.lambda}{the data set for the extrapolation,}
#' \item{variance.asymptotic}{the asymptotic variance estimates,}
#' \item{theta}{all estimated coefficients for each lambda and B,}
#' ...
#'
#' @details
#' If \code{mc.matrix} is a function the first argument of that function must
#' be the whole dataset used in the naive model, the second argument must be
#' the exponent (lambda) for the misclassification. The function must return
#' a \code{data.frame} containing the misclassified \code{SIMEXvariable}. An
#' example can be found below.
#'
#' Asymptotic variance estimation is only implemented for \code{lm} and \code{glm}
#'
#' The loglinear fit has the form \code{g(lambda, GAMMA) = exp(gamma0 + gamma1 * lambda)}.
#' It is realized via the \code{log()} function. To avoid negative values the
#' minimum +1 of the dataset is added and after the prediction later substracted
#' \code{exp(predict(...)) - min(data) - 1}.
#'
#' The 'log2' fit is fitted via the \code{nls()} function for direct fitting of
#' the model \code{y ~ exp(gamma.0 + gamma.1 * lambda)}. As starting values the
#' results of a LS-fit to a linear model with a log transformed response are used.
#' If \code{nls} does not converge, the model with the starting values is returned.
#'
#' \code{refit()} refits the object with a different extrapolation function.
#'
#' @references
#'Küchenhoff, H., Mwalili, S. M.  and Lesaffre, E. (2006) A general method for
#'dealing with misclassification in regression: The Misclassification SIMEX.
#'\emph{Biometrics}, \bold{62}, 85 -- 96
#'
#' Küchenhoff, H., Lederer, W. and E. Lesaffre. (2006) Asymptotic Variance Estimation
#' for the Misclassification SIMEX.
#' \emph{Computational Statistics and Data Analysis}, \bold{51}, 6197 -- 6211
#'
#' Lederer, W. and Küchenhoff, H. (2006) A short introduction to the SIMEX and MCSIMEX. \emph{R News}, \bold{6(4)}, 26--31
#'
#' @author Wolfgang Lederer, \email{wolfgang.lederer@gmail.com}
#' @seealso \code{\link[simex]{misclass}}, \code{\link[simex]{simex}}
#'
#' @examples
#' x <- rnorm(200, 0, 1.142)
#' z <- rnorm(200, 0, 2)
#' y <- factor(rbinom(200, 1, (1 / (1 + exp(-1 * (-2 + 1.5 * x -0.5 * z))))))
#' Pi <- matrix(data = c(0.9, 0.1, 0.3, 0.7), nrow = 2, byrow = FALSE)
#' dimnames(Pi) <- list(levels(y), levels(y))
#' ystar <- misclass(data.frame(y), list(y = Pi), k = 1)[, 1]
#' naive.model <- glm(ystar ~ x + z, family = binomial, x = TRUE, y = TRUE)
#' true.model  <- glm(y ~ x + z, family = binomial)
#' simex.model <- mcsimex(naive.model, mc.matrix = Pi, SIMEXvariable = "ystar")
#'
#' op <- par(mfrow = c(2, 3))
#' invisible(lapply(simex.model$theta, boxplot, notch = TRUE, outline = FALSE,
#'                  names = c(0.5, 1, 1.5, 2)))
#'                  plot(simex.model)
#' simex.model2 <- refit(simex.model, "line")
#' plot(simex.model2)
#' par(op)
#'
#' # example for a function which can be supplied to the function mcsimex()
#' # "ystar" is the variable which is to be misclassified
#' # using the example above
#' \dontrun{
#' my.misclass <- function (datas, k) {
#'     ystar <- datas$"ystar"
#'     p1 <- matrix(data = c(0.75, 0.25, 0.25, 0.75), nrow = 2, byrow = FALSE)
#'     colnames(p1) <- levels(ystar)
#'     rownames(p1) <- levels(ystar)
#'     p0 <- matrix(data = c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = FALSE)
#'
#'     colnames(p0) <- levels(ystar)
#'     rownames(p0) <- levels(ystar)
#'     ystar[datas$x < 0] <-
#'     misclass(data.frame(ystar = ystar[datas$x < 0]), list(ystar = p1), k = k)[, 1]
#'     ystar[datas$x > 0] <-
#'     misclass(data.frame(ystar = ystar[datas$x > 0]), list(ystar = p0), k = k)[, 1]
#'     ystar <- factor(ystar)
#'     return(data.frame(ystar))}
#'
#' simex.model.differential <- mcsimex(naive.model, mc.matrix = "my.misclass", SIMEXvariable = "ystar")
#' }
#'
#' @keywords models
#' @export



mcsimex <- function(model,
                    SIMEXvariable,
                    mc.matrix,
                    lambda = c(0.5, 1,1.5, 2),
                    B = 100,
                    fitting.method = "quadratic",
                    jackknife.estimation = "quadratic",
                    asymptotic = TRUE) {
  fitting.method <- substr(fitting.method, 1, 4)
  if (!any(fitting.method == c("quad", "line", "nonl", "logl", "log2"))) {
    warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
    fitting.method <- "quad"
  }
  if (jackknife.estimation != FALSE)
    jackknife.estimation <- substr(jackknife.estimation, 1, 4)
  if (!any(jackknife.estimation == c("quad", "line", "nonl", "logl", FALSE))) {
    warning("Fitting Method (jackknife) not implemented. Using: quadratic",
            call. = FALSE)
    jackknife.estimation <- "quad"
  }
  if (any(lambda <= 0)) {
    warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
            call. = FALSE)
    lambda <- lambda[lambda >= 0]
  }
  if (class(model)[1] == "polr" && !any(names(model)) == "Hessian")
    stop("The option Hessian must be enabled in the naive model", call. = FALSE)
  if (!any(names(model) == "x") && asymptotic && class(model)[1] != "polr")
    stop("The option x must be enabled in the naive model for asymptotic variance estimation",
         call. = FALSE)
  if (is.matrix(mc.matrix)) {
    mc.matrix <- list(mc.matrix)
    names(mc.matrix) <- SIMEXvariable
  }
  if (is.list(mc.matrix)) {
    if (!all(check.mc.matrix(mc.matrix)))
      stop("mc.matrix may contain negative values for exponents smaller than 1")
    if (length(mc.matrix) != length(SIMEXvariable))
      stop("mc.matrix and SIMEXvariable do not match")
  }
  if (any(!sapply(as.data.frame(model$model), is.factor)[SIMEXvariable]))
    stop("SIMEXvariable must be a factor")
  cl <- match.call()
  ncoef <- length(model$coefficients)
  if (class(model)[1] == "polr")
    ncoef <- ncoef + length(model$zeta)
  ndes <- length(model$y)
  nlambda <- length(lambda)
  if (class(model)[1] == "polr")
    p.names <- c(names(coef(model)), names(model$zeta))
  else
    p.names <- names(coef(model))
  factors <- lapply(SIMEXvariable, levels)
  estimates <- matrix(data = NA, length(lambda) + 1, length(model$coefficients))
  theta <- matrix(data = NA, B, ncoef)
  colnames(theta) <- p.names
  theta.all <- vector(mode = "list", nlambda)
  if (jackknife.estimation != FALSE) {
    var.exp <- list()
    var.exp[[1]] <- extract.covmat(model)
  }
  if (asymptotic) {
    psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
    psi <- residuals(model, type = "response") * model$x
    PSI <- psi
    am <- list()
    xi <- model$x
    a <- list()
    dh <- model$family$mu.eta(model$linear.predictors)
    for (k in 1:ndes) a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k, ])
    a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
    ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
    am[[1]] <- -ab/ndes
  }
  if (class(model)[1] == "polr")
    estimates[1, ] <- c(model$coefficients, model$zeta)
  else
    estimates[1, ] <- model$coefficients
  for (i in 1:length(lambda)) {
    if (jackknife.estimation != FALSE)
      variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
    if (asymptotic) {
      psi <- matrix(0, ncol = ncoef, nrow = ndes)
      a <- list()
      for (k in 1:ndes) a[[k]] <- matrix(0, nrow = ncoef, ncol = ncoef)
    }
    for (j in 1:B) {
      SIMEXdata <- data.frame(model$model)
      # doing the misclassification
      SIMEXv <- data.frame(SIMEXdata[, SIMEXvariable])
      colnames(SIMEXv) <- SIMEXvariable
      if (is.character(mc.matrix)) {
        SIMEXdata[, SIMEXvariable] <- eval(call(mc.matrix, SIMEXdata, lambda[i]))
      } else {
        SIMEXdata[, SIMEXvariable] <- misclass(SIMEXv, mc.matrix, lambda[i])
      }
      # updating the model and calculating the estimates
      model.SIMEX <- update(model, data = data.frame(SIMEXdata))
      if (class(model)[1] == "polr")
        theta[j, ] <- c(model.SIMEX$coefficients, model.SIMEX$zeta)
      else
        theta[j, ] <- model.SIMEX$coefficients
      if (jackknife.estimation != FALSE) {
        variance.est <- variance.est + extract.covmat(model.SIMEX)
      }
      if (asymptotic) {
        xi <- model.SIMEX$x
        psi <- psi + (residuals(model.SIMEX, type = "response") * xi)
        dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
        for (k in 1:ndes) a[[k]] <- a[[k]] - dh[k] * xi[k, ] %*% t(xi[k, ])
      }
    }
    # taking the mean of the estimate -> SIMEX estimate
    estimates[i + 1, ] <- colMeans(theta)
    theta.all[[i]] <- theta
    if (jackknife.estimation != FALSE) {
      variance.est <- variance.est/B
      s2 <- cov(theta)
      var.exp[[i + 1]] <- variance.est - s2
    }
    if (asymptotic) {
      xiB <- psi/B
      PSI <- cbind(PSI, xiB)
      a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
      ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
      am[[i + 1]] <- ab/(B * ndes)
    }
  }
  lambda <- c(0, lambda)
  colnames(estimates) <- p.names
  # fitting the extrapolation function
  switch(fitting.method,
         quad = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
         line = extrapolation <- lm(estimates ~ lambda),
         logl = extrapolation <- lm(I(log(t(t(estimates) + (abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0)))) ~ lambda),
         log2 = extrapolation <- fit.logl(lambda, p.names, estimates),
         nonl = extrapolation <- fit.nls(lambda, p.names, estimates))
  if (any(class(extrapolation) == "lm") && fitting.method == "log2")
    fitting.method <- "logl"
  # predicting the SIMEX estimate
  SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
  switch(fitting.method, quad = SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1)),
         line = SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1)),
         nonl = for (i in 1:length(p.names)) SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]], newdata = data.frame(lambda = -1)),
         logl = SIMEX.estimate <- exp(predict(extrapolation, newdata = data.frame(lambda = -1))) - (abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0),
         log2 = for (i in 1:length(p.names)) SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]], newdata = data.frame(lambda = -1)) - ((abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0))[i])

  if (jackknife.estimation != FALSE) {
    variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2, byrow = TRUE)
    switch(jackknife.estimation,
           quad = extrapolation.variance <- lm(variance.jackknife ~ lambda + I(lambda^2)),
           line = extrapolation.variance <- lm(variance.jackknife ~ lambda),
           logl = extrapolation.variance <- lm(I(log(t(t(variance.jackknife) + (abs(apply(variance.jackknife, 2, min)) + 1) * (apply(variance.jackknife, 2, min) <= 0)))) ~ lambda),
           nonl = extrapolation.variance <- fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife))
    variance.jackknife2 <- vector("numeric", ncoef^2)
    switch(jackknife.estimation,
           nonl = for (i in 1:NCOL(variance.jackknife)) variance.jackknife2[i] <- predict(extrapolation.variance[[i]], newdata = data.frame(lambda = -1)),
           quad = variance.jackknife2 <- predict(extrapolation.variance, newdata = data.frame(lambda = -1)),
           line = variance.jackknife2 <- predict(extrapolation.variance, newdata = data.frame(lambda = -1)),
           logl = variance.jackknife2 <- exp(predict(extrapolation.variance, newdata = data.frame(lambda = -1))) - (abs(apply(variance.jackknife, 2, min)) + 1) * (apply(variance.jackknife, 2, min) <= 0))
    variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
    variance.jackknife.lambda <- cbind(c(-1, lambda), variance.jackknife)
    variance.jackknife <- matrix(variance.jackknife[1, ], nrow = ncoef,
                                 ncol = ncoef, byrow = TRUE)
    dimnames(variance.jackknife) <- list(p.names, p.names)
  }
  if (asymptotic) {
    c11 <- cov(PSI)
    a11 <- diag.block(am)
    a11.inv <- solve(a11)
    sigma <- a11.inv %*% c11 %*% t(a11.inv)
    s <- construct.s(ncoef, lambda, fitting.method, extrapolation)
    d.inv <- solve(s %*% t(s))
    sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
    g <- list()
    switch(fitting.method,
           quad = g <- c(1, -1, 1), line = g <- c(1, -1),
           logl = for (i in 1:ncoef) g[[i]] <- c(exp(coef(extrapolation)[1,i] - coef(extrapolation)[2, i]), -exp(coef(extrapolation)[1,i] - coef(extrapolation)[2, i])),
           log2 = for (i in 1:ncoef) g[[i]] <- c(exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2]), -exp(coef(extrapolation[[i]])[1] - coef(extrapolation[[i]])[2])), nonl = for (i in 1:ncoef) g[[i]] <- c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1, coef(extrapolation[[i]])[2]/(coef(extrapolation[[i]])[3] - 1)^2))
    g <- diag.block(g, ncoef)
    variance.asymptotic <- (t(g) %*% sigma.gamma %*% g)/ndes
    dimnames(variance.asymptotic) <- list(p.names, p.names)
  }
  # creating class 'mcsimex'
  theta <- matrix(unlist(theta.all), nrow = B)
  theta.all <- list()
  for (i in 1:ncoef) theta.all[[p.names[i]]] <- data.frame(theta[, seq(i, ncoef * nlambda, by = ncoef)])
  z <- cbind(lambda, estimates)
  z <- rbind(c(-1, SIMEX.estimate), z)  # returning the estimated values
  colnames(z) <- c("lambda", p.names)
  if(class(model)[1] == "polr")
    coefs <- z[1, -1][1:length(coef(model))]
  else
    coefs <- z[1, -1]
  erg <- list(coefficients = coefs, SIMEX.estimates = z, lambda = lambda,
              model = model, mc.matrix = mc.matrix, B = B, extrapolation = extrapolation,
              fitting.method = fitting.method, SIMEXvariable = SIMEXvariable,
              call = cl, theta = theta.all)
  class(erg) <- ("mcsimex")
  if (class(model)[1] == "polr")
    erg$zeta <- z[1,-1][(length(coef(model))+1):length(z[1,-1])]
  type <- 'response'
  if (class(model)[1] == "polr")
    type <- 'probs'
  fitted.values <- predict(erg, newdata = model$model[, -1, drop = FALSE],
                           type = type)
  erg$fitted.values <- fitted.values
  if (is.factor(model$model[, 1]))
    erg$residuals <- as.numeric(levels(model$model[, 1]))[model$model[, 1]] - fitted.values
  else erg$residuals <- model$model[, 1] - fitted.values

  if (jackknife.estimation != FALSE) {
    erg$extrapolation.variance <- extrapolation.variance
    erg$variance.jackknife <- variance.jackknife
    erg$variance.jackknife.lambda <- variance.jackknife.lambda
  }
  if (asymptotic) {
    erg$PSI <- PSI
    erg$c11 <- c11
    erg$a11 <- a11
    erg$sigma <- sigma
    erg$sigma.gamma <- sigma.gamma
    erg$g <- g
    erg$s <- s
    erg$variance.asymptotic <- variance.asymptotic
  }
  return(erg)
}

#' @describeIn mcsimex Plots of the simulation and extrapolation
#' @export

plot.mcsimex <- function(x,
                         xlab = expression((1 + lambda)),
                         ylab = colnames(b[, -1]),
                         ask = FALSE,
                         show = rep(TRUE, NCOL(b) - 1),
                         ...) {
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(...)
  if (ask)
    par(ask = TRUE)
  p.names <- names(coef(x))
  b <- x$SIMEX.estimates
  a <- seq(-1, max(b[, 1]), by = 0.01)
  d <- matrix(data = NA, nrow = length(a), ncol = NCOL(b) - 1)
  switch(x$fitting.method,
         quad = d <- predict(x$extrapolation, newdata = data.frame(lambda = a)),
         line = d <- predict(x$extrapolation, newdata = data.frame(lambda = a)),
         nonl = for (i in 1:length(p.names)) d[, i] <- predict(x$extrapolation[[p.names[i]]],
                                                               newdata = data.frame(lambda = a)),
         log2 = for (i in 1:length(p.names)) d[, i] <- predict(x$extrapolation[[p.names[i]]], newdata = data.frame(lambda = a)) -
                 ((abs(apply(x$SIMEX.estimates[-1, -1], 2, min)) + 1) * (apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0))[i],
         logl = d <- t(t(exp(predict(x$extrapolation, newdata = data.frame(lambda = a)))) - ((abs(apply(x$SIMEX.estimates[-1,
                  -1], 2, min)) + 1) * (apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0))))
  for (i in 2:NCOL(b)) {
    if (show[i - 1]) {
      plot(b[, 1] + 1, b[, i], xlab = xlab, ylab = ylab[i - 1], type = "n")
      points(b[-1, 1] + 1, b[-1, i], pch = 19)
      points(b[1, 1] + 1, b[1, i])
      lines(a[a > 0] + 1, d[a > 0, (i - 1)])
      lines(a[a < 0] + 1, d[a < 0, (i - 1)], lty = 2)
    }
  }
}

#' @describeIn mcsimex predict with mcsimex correction
#' @export
predict.mcsimex <- function(object, newdata, ...) {
  new.object <- object$model
  new.object$coefficients <- object$coefficients
  if (class(new.object)[1] == "polr") {
    new.object$zeta <- object$zeta
    new.object$fitted.values <- object$fitted.values
  }
  if (missing(newdata)) {
    predict(new.object, ...)
  } else {
    predict(new.object, newdata = data.frame(newdata), ...)
  }
}

#' @describeIn mcsimex Nice printing
#' @export
print.mcsimex <- function(x, digits = max(3, getOption("digits") - 3),
                          ...) {
  cat("\nNaive model:\n", deparse(x$model$call), "\n", sep = "")
  cat("\nSIMEX-Variables: ")
  cat(x$SIMEXvariable, sep = ", ")
  cat("\nNumber of Simulations: ", paste(x$B), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2,
                  quote = FALSE)
  } else cat("No coefficients\n")
  if (length(x$zeta)) {
    cat("Intercepts:\n")
    print.default(format(x$zeta, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  cat("\n")
  return(invisible(x))
}

#' @describeIn mcsimex Print summary nicely
#' @export
print.summary.mcsimex <- function(x, digits = max(3, getOption("digits") -
                                                    3), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNaive model: \n")
  print(x$naive.model)
  cat("\nSimex variable :", x$SIMEXvariable, "\n")
  cat("Misclassification matrix: \n")
  if (is.character(x$mc.matrix))
    print(x$mc.matrix) else lapply(x$mc.matrix, print)
  cat("\nNumber of iterations: ", x$B, "\n")
  cat("\nResiduals: \n")
  print(summary(x$residuals), digits)
  cat("\nCoefficients: \n")
  if (any(names(x$coefficients) == "asymptotic")) {
    cat("\nAsymptotic variance: \n")
    printCoefmat(x$coefficients$asymptotic, digits = digits)
  }
  if (any(names(x$coefficients) == "jackknife")) {
    cat("\nJackknife variance: \n")
    printCoefmat(x$coefficients$jackknife, digits = digits)
  }
  return(invisible(x))
}

#' @describeIn mcsimex Summary for mcsimex
#' @export
summary.mcsimex <- function(object, ...) {
  if (class(object$model)[1] == "polr")
    est <- c(coef(object), object$zeta)
  else
    est <- coef(object)
  p.names <- names(est)
  est.table <- list()
  n <- length(resid(object))
  p <- length(p.names)
  rdf <- n - p
  if (any(names(object) == "variance.jackknife")) {
    se <- sqrt(diag(object$variance.jackknife))
    tval <- est/se
    pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
    est.table[["jackknife"]] <- cbind(est, se, tval, pval)
    dimnames(est.table[["jackknife"]]) <- list(p.names, c("Estimate",
                                                          "Std. Error", "t value", "Pr(>|t|)"))
  }
  if (any(names(object) == "variance.asymptotic")) {
    se <- sqrt(diag(object$variance.asymptotic))
    tval <- est/se
    pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
    est.table[["asymptotic"]] <- cbind(est, se, tval, pval)
    dimnames(est.table[["asymptotic"]]) <- list(p.names, c("Estimate",
                                                           "Std. Error", "t value", "Pr(>|t|)"))
  }
  ans <- list()
  class(ans) <- "summary.mcsimex"
  ans$coefficients <- est.table
  ans$residuals <- resid(object)
  ans$call <- object$call
  ans$B <- object$B
  ans$naive.model <- object$model$call
  ans$SIMEXvariable <- object$SIMEXvariable
  ans$mc.matrix <- object$mc.matrix
  return(ans)
}

#' @describeIn mcsimex Refits the model with a different extrapolation function
#' @export
refit.mcsimex <- function(object, fitting.method = "quadratic", jackknife.estimation = "quadratic",
                          asymptotic = TRUE, ...)
              .refit(object, fitting.method = fitting.method, jackknife.estimation = jackknife.estimation,
                     asymptotic = asymptotic, ...)





