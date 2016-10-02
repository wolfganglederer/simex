#' Measurement error in models using SIMEX
#'
#' Implementation of the SIMEX algorithm for measurement error models according to Cook and Stefanski
#'
#' @aliases simex print.simex summary.simex print.summary.simex plot.simex predict.simex refit refit.simex
#'
#' @param model the naive model
#' @param SIMEXvariable character or vector of characters containing the names of the variables with measurement error
#' @param measurement.error given standard deviations of measurement errors. In
#' case of homoskedastic measurement error it is a matrix with dimension
#' 1x\code{length(SIMEXvariable)}. In case of heteroskedastic error for at least one
#' SIMEXvariable it is a matrix of dimension nx
#' @param lambda vector of lambdas for which the simulation step should be done (without 0)
#' @param B number of iterations for each lambda
#' @param fitting.method fitting method for the extrapolation. \code{linear}, \code{quadratic},
#' \code{nonlinear} are implemented. (first 4 letters are enough)
#' @param jackknife.estimation specifying the extrapolation method for jackknife
#' variance estimation. Can be set to \code{FALSE} if it should not be performed
#' @param asymptotic logical, indicating if asymptotic variance estimation should
#' be done, in the naive model the option \code{x = TRUE} has to be set
#' @param x object of class 'simex'
#' @param digits number of digits to be printed
#' @param object of class 'simex'
#' @param xlab optional name for the X-Axis
#' @param ylab vector containing the names for the Y-Axis
#' @param ask logical. If \code{TRUE}, the user is asked for input, before a new figure is drawn
#' @param show vector of logicals indicating for wich variables a plot should be produced
#' @param newdata optionally, a data frame in which to look for
#' variables with which to predict. If omitted, the fitted linear predictors are used
#' @param \dots arguments passed to other functions
#'
#' @details
#'
#' Nonlinear is implemented as described in Cook and Stefanski, but is numerically
#' instable. It is not advisable to use this feature. If a nonlinear extrapolation
#' is desired please use the \code{refit()} method.
#'
#' Asymptotic is only implemented for naive models of class \code{lm} or \code{glm} with homoscedastic measurement error.
#'
#' \code{refit()} refits the object with a different extrapolation function.
#'
#' @return
#' An object of class 'simex' which contains:
#'
#' \item{coefficients}{the corrected coefficients of the SIMEX model,}
#' \item{SIMEX.estimates}{the estimates for every lambda,}
#' \item{model}{the naive model,}
#' \item{measurement.error}{the known error standard deviations,}
#' \item{B}{the number of iterations,}
#' \item{extrapolation}{the model object of the extrapolation step,}
#' \item{fitting.method}{the fitting method used in the extrapolation step,}
#' \item{residuals}{the residuals of the main model,}
#' \item{fitted.values}{the fitted values of the main model,}
#' \item{call}{the function call,}
#' \item{variance.jackknife}{the jackknife variance estimate,}
#' \item{extrapolation.variance}{the model object of the variance extrapolation,}
#' \item{variance.jackknife.lambda}{the data set for the extrapolation,}
#' \item{variance.asymptotic}{the asymptotic variance estimates,}
#' \item{theta}{the estimates for every B and lambda,}
#' ...
#'
#' @references
#' Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in
#' parametric measurement error models. \emph{Journal of the American Statistical
#' Association}, \bold{89}, 1314 -- 1328
#'
#' Carroll, R.J., Küchenhoff, H., Lombard, F. and Stefanski L.A. (1996) Asymptotics
#' for the SIMEX estimator in nonlinear measurement error models. \emph{Journal
#'   of the American Statistical Association}, \bold{91}, 242 -- 250
#'
#' Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006).
#' \emph{Measurement error in nonlinear models: A modern perspective.}, Second
#' Edition. London: Chapman and Hall.
#'
#' Lederer, W. and Küchenhoff, H. (2006) A short introduction to the SIMEX and MCSIMEX. \emph{R News}, \bold{6(4)}, 26--31
#'
#' @author Wolfgang Lederer,\email{wolfgang.lederer@gmail.com}
#' @author Heidi Seibold,\email{heidi.bold@gmail.com}
#'
#' @seealso \code{\link[simex]{mcsimex}} for discrete data with misclassification,
#'  \code{\link[stats]{lm}}, \code{\link[stats]{glm}}
#'
#' @examples
#' ## Seed
#' set.seed(49494)
#'
#' ## simulating the measurement error standard deviations
#' sd_me <- 0.3
#' sd_me2 <- 0.4
#' temp <- runif(100, min = 0, max = 0.6)
#' sd_me_het1 <- sort(temp)
#' temp2 <- rnorm(100, sd = 0.1)
#' sd_me_het2 <- abs(sd_me_het1 + temp2)
#'
#' ## simulating the independent variables x (real and with measurement error):
#'
#' x_real <- rnorm(100)
#' x_real2 <- rpois(100, lambda = 2)
#' x_real3 <- -4*x_real + runif(100, min = -10, max = 10)  # correlated to x_real
#'
#' x_measured <- x_real + sd_me * rnorm(100)
#' x_measured2 <- x_real2 + sd_me2 * rnorm(100)
#' x_het1 <- x_real + sd_me_het1 * rnorm(100)
#' x_het2 <- x_real3 + sd_me_het2 * rnorm(100)
#'
#' ## calculating dependent variable y:
#' y <- x_real + rnorm(100, sd = 0.05)
#' y2 <- x_real + 2*x_real2 + rnorm(100, sd = 0.08)
#' y3 <- x_real + 2*x_real3 + rnorm(100, sd = 0.08)
#'
#' ### one variable with homoscedastic measurement error
#' (model_real <- lm(y ~ x_real))
#'
#' (model_naiv <- lm(y ~ x_measured, x = TRUE))
#'
#' (model_simex <- simex(model_naiv, SIMEXvariable = "x_measured", measurement.error = sd_me))
#' plot(model_simex)
#'
#' ### two variables with homoscedastic measurement errors
#' (model_real2 <- lm(y2 ~ x_real + x_real2))
#' (model_naiv2 <- lm(y2 ~ x_measured + x_measured2, x = TRUE))
#' (model_simex2 <- simex(model_naiv2, SIMEXvariable = c("x_measured", "x_measured2"),
#'          measurement.error = cbind(sd_me, sd_me2)))
#'
#' plot(model_simex2)
#'
#' \dontrun{
#' ### one variable with increasing heteroscedastic measurement error
#' model_real
#'
#' (mod_naiv1 <- lm(y ~ x_het1, x = TRUE))
#' (mod_simex1 <- simex(mod_naiv1, SIMEXvariable = "x_het1",
#'                 measurement.error = sd_me_het1, asymptotic = FALSE))
#'
#' plot(mod_simex1)
#'
#' ### two correlated variables with heteroscedastic measurement errors
#' (model_real3 <- lm(y3 ~ x_real + x_real3))
#' (mod_naiv2 <- lm(y3 ~ x_het1 + x_het2, x = TRUE))
#' (mod_simex2 <- simex(mod_naiv2, SIMEXvariable = c("x_het1", "x_het2"),
#'               measurement.error = cbind(sd_me_het1, sd_me_het2), asymptotic = FALSE))
#'
#' plot(mod_simex2)
#'
#' ### two variables, one with homoscedastic, one with heteroscedastic measurement error
#' model_real2
#' (mod_naiv3 <- lm(y2 ~ x_measured + x_het2, x = TRUE))
#' (mod_simex3 <- simex(mod_naiv3, SIMEXvariable = c("x_measured", "x_het2"),
#'                     measurement.error = cbind(sd_me, sd_me_het2), asymptotic = FALSE))
#'
#' ### glm: two variables, one with homoscedastic, one with heteroscedastic measurement error
#' t <- x_real + 2*x_real2 + rnorm(100, sd = 0.01)
#' g <- 1 / (1 + exp(t))
#' u <- runif(100)
#' ybin <- as.numeric(u < g)
#'
#' (logit_real <- glm(ybin ~ x_real + x_real2, family = binomial))
#' (logit_naiv <- glm(ybin ~ x_measured + x_het2, x = TRUE, family = binomial))
#' (logit_simex <- simex(logit_naiv, SIMEXvariable = c("x_measured", "x_het2"),
#'                     measurement.error = cbind(sd_me, sd_me_het2), asymptotic = FALSE))
#'
#' summary(logit_simex)
#' print(logit_simex)
#' plot(logit_simex)
#'
#' ### polr: two variables, one with homoscedastic, one with heteroscedastic measurement error
#'
#' if(require("MASS")) {# Requires MASS
#' yerr <- jitter(y, amount=1)
#' yfactor <- cut(yerr, 3, ordered_result=TRUE)
#'
#' (polr_real <- polr(yfactor ~ x_real + x_real2))
#' (polr_naiv <- polr(yfactor ~ x_measured + x_het2, Hess = TRUE))
#' (polr_simex <- simex(polr_naiv, SIMEXvariable = c("x_measured", "x_het2"),
#'                     measurement.error = cbind(sd_me, sd_me_het2), asymptotic = FALSE))
#'
#' summary(polr_simex)
#' print(polr_simex)
#' plot(polr_simex)
#' }
#' }
#'
#' @keywords models
#'
#' @export

simex <-
  function (
    model,                          # the naive model
    SIMEXvariable,                  # vector of names for the SIMEXvariables
    measurement.error,              # vector of standard deviations of the measurement errors
    lambda = c(0.5, 1, 1.5, 2),     # vector of lambda which contains the values for the ...
    B = 100,                        # numeric: number of simulations to be made in each step
    fitting.method = "quadratic",   # fitting method for the extrapolation step
    jackknife.estimation = "quadratic", # extrapolation function for the variance estimation
    asymptotic = TRUE)              # logical: do asymptotic variance estimation
  {
    fitting.method <- substr(fitting.method, 1, 4)
    if (!any(fitting.method == c("quad", "line", "nonl"))) {
      warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
      fitting.method <- "quad"
    }
    if (jackknife.estimation != FALSE)
      jackknife.estimation <- substr(jackknife.estimation, 1, 4)
    if (!any(jackknife.estimation == c("quad", "line", "nonl", FALSE))) {
      warning("Fitting Method (jackknife) not implemented. Using: quadratic",
              call. = FALSE)
      jackknife.estimation <- "quad"
    }
    if (!is.character(SIMEXvariable))
      stop("SIMEXvariable must be character", call. = FALSE)
    if (any(lambda <= 0)) {
      warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
              call. = FALSE)
      lambda <- lambda[lambda >= 0]
    }
    if (class(model)[1] == "polr" && !any(names(model) == "Hessian"))
      stop("The option Hessian must be enabled in the naive model", call. = FALSE)
    if (!any(names(model) == "x") && asymptotic && class(model)[1] != "polr")
      stop("The option x must be enabled in the naive model for asymptotic variance estimation",
           call. = FALSE)
    #**Heidi**#
    measurement.error <- as.matrix(measurement.error)
    SIMEXdata <- model$model #das war vorher nach "for (j in 1:B)"
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) == 1){
      measurement.error <- matrix(measurement.error, nrow = NROW(SIMEXdata), ncol = NCOL(measurement.error), byrow = TRUE)
    }
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) != 1){
      stop("NROW(measurement.error) must be either 1 or must take the number of rows of the data used.",
           call. = FALSE)
    }
    # may now be 0
    if (any(measurement.error < 0))
      stop("measurement.error is negative", call. = FALSE)
    # but not constant 0
    any0 <- ( apply(measurement.error, 2, all.equal, current = rep(0, times = NROW(measurement.error))) == TRUE )
    if ( sum(any0) > 0 )
      stop("measurement.error is constant 0 in column(s) ", which(any0), call. = FALSE)

    # Abbruch falls das Modell nicht lm oder glm bzw. heterosk. me und asymptotic = TRUE
    if ( asymptotic == TRUE & ((class(model)[1] != "glm" & class(model)[1] != "lm") | dim(unique(measurement.error))[1] != 1) )
      stop("Asymptotic is only implemented for naive models of class lm or glm with homoscedastic measurement error.")
    cl <- match.call()
    # defining the vector for the solutions of the simulations
    ncoef <- length(model$coefficients)
    if (class(model)[1] == "polr")
      ncoef <- ncoef + length(model$zeta)
    ndes <- dim(model$model)[1]
    if (class(model)[1] == "polr")
      p.names <- c(names(coef(model)), names(model$zeta))
    else
      p.names <- names(coef(model))
    nlambda <- length(lambda)
    estimates <- matrix(data = NA, nlambda + 1, ncoef) # +1 because "0" will be added
    theta <- matrix(data = NA, B, ncoef)
    colnames(theta)<- p.names
    theta.all <-  vector(mode = "list", nlambda)
    if (jackknife.estimation != FALSE) {
      var.exp <- list()
      var.exp[[1]] <- extract.covmat(model)
    }
    if (asymptotic) {
      psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
      psi <- resid(model, type = "response") * model$x
      PSI <- psi
      am <- list()
      a <- list()
      xi <- model$x
      dh <- rep(1, ndes)
      if (class(model)[1] == "glm")
        dh <- model$family$mu.eta(model$linear.predictors)
      for (k in 1:ndes)
        a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k, ])
      a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
      ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
      am[[1]] <- -ab / ndes
      a <- list()
    }
    # assigning the naive estimator
    if (class(model)[1] == "polr")
      estimates[1, ] <- c(model$coefficients, model$zeta)
    else
      estimates[1, ] <- model$coefficients
    # The Simulation step
    # outer loop doing the simulations for each lambda
    for (i in 1:length(lambda)) {
      if (jackknife.estimation != FALSE)
        variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
      if (asymptotic) {
        psi <- matrix(0, ncol = ncoef, nrow = ndes)
        a <- list()
        for(k in 1:ndes)
          a[[k]] <- matrix(0, nrow = ncoef, ncol = ncoef)
      }
      # inner loop, doing the simulations
      for (j in 1:B) {
        SIMEXdata <- model$model
        epsilon <- matrix(rnorm(n=NROW(SIMEXdata)*length(SIMEXvariable))
                          ,ncol= length(SIMEXvariable)
                          ,nrow =  NROW(SIMEXdata))
        # and adding the random error
        SIMEXdata[, SIMEXvariable] <- SIMEXdata[, SIMEXvariable] +
          ( sqrt(lambda[i]) *  epsilon * measurement.error  )
        # updating the model and calculating the estimate
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
          psi <- psi + (resid(model.SIMEX, type = "response") * xi)
          dh <- rep(1, ndes)
          if (class(model)[1] == "glm")
            dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
          for (k in 1:ndes)
            a[[k]] <- a[[k]] - dh[k] * xi[k, ] %*% t(xi[k, ])
        }
      }
      # taking the mean of the estimate -> SIMEX estimate
      estimates[i + 1, ] <- colMeans(theta)
      theta.all[[i]] <- theta
      # Variance estimation via the Jackknife
      if (jackknife.estimation != FALSE) {
        variance.est <- variance.est / B
        s2 <- cov(theta)
        var.exp[[i + 1]]<- variance.est - s2
      }
      if (asymptotic) {
        xiB <- psi / B
        PSI <- cbind(PSI, xiB)
        a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
        ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
        am[[i + 1]] <- ab / (B * ndes)
      }
    }
    # extrapolation step
    SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
    colnames(estimates) <- p.names
    lambda <- c(0, lambda)
    # fitting the extrapolation function
    switch(fitting.method,
           "quad" = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
           "line" = extrapolation <- lm(estimates ~ lambda),
           "nonl" = extrapolation <- fit.nls(lambda, p.names, estimates)
    )
    # predicting the SIMEX estimate
    if (fitting.method == "nonl") {
      for (i in 1:length(p.names))
        SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
                                     newdata = data.frame(lambda = -1))
    } else {
      SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1))
    }
    ######  Jackknife Estimation
    if (jackknife.estimation != FALSE) {
      variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2 , byrow = TRUE)
      switch(jackknife.estimation,
             "quad" = extrapolation.variance <-
               lm(variance.jackknife ~ lambda + I(lambda^2)),
             "line" = extrapolation.variance <- lm(variance.jackknife ~ lambda),
             "nonl" = extrapolation.variance <-
               fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife)
      )
      variance.jackknife2 <- vector("numeric", ncoef^2)
      switch(jackknife.estimation,
             "nonl"= for(i in 1:NCOL(variance.jackknife))
               variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
                                                 newdata = data.frame(lambda = -1)),
             "quad"= variance.jackknife2 <- predict(extrapolation.variance,
                                                    newdata = data.frame(lambda = -1)),
             "line"= variance.jackknife2 <- predict(extrapolation.variance,
                                                    newdata = data.frame(lambda = -1))
      )
      variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
      variance.jackknife.lambda <- cbind(c(-1,lambda), variance.jackknife)
      variance.jackknife <- matrix(variance.jackknife[1, ], nrow = ncoef,
                                   ncol = ncoef, byrow = TRUE)
      dimnames(variance.jackknife) <- list(p.names, p.names)
    }
    # Asymptotic estimation
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
             "quad" = g <- c(1, -1, 1),
             "line" = g <- c(1, -1),
             "nonl" = for(i in 1:ncoef)
               g[[i]] <- c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1,
                           coef(extrapolation[[i]])[2] / (coef(extrapolation[[i]])[3] - 1)^2)
      )
      g <- diag.block(g, ncoef)
      variance.asymptotic <- (t(g) %*% sigma.gamma %*% g) / ndes
      dimnames(variance.asymptotic) <- list(p.names, p.names)
    }
    # creating class "simex"
    theta <- matrix(unlist(theta.all), nrow = B)
    theta.all <- list()
    for (i in 1:ncoef)
      theta.all[[p.names[i]]] <-
      data.frame(theta[, seq(i, ncoef * nlambda, by = ncoef)])
    z <- cbind(lambda, estimates)
    z <- rbind(c(-1, SIMEX.estimate), z) # returning the estimated values
    colnames(z) <- c("lambda", p.names)
    if(class(model)[1] == "polr")
      coefs <- z[1, -1][1:length(coef(model))]
    else
      coefs <- z[1, -1]
    erg <- list(
      coefficients = coefs, # SIMEX corrected coefficients
      SIMEX.estimates = z, # all thetas as a matrix
      lambda = lambda, # vector for the values for lambda
      model = model, # the naive model
      measurement.error = measurement.error, # vector of values of measurement.error
      B = B, # number of Simulations
      extrapolation = extrapolation, # model of the extrapolation
      fitting.method = fitting.method,# which fitting method was used
      SIMEXvariable = SIMEXvariable,
      theta = theta.all,
      call = cl
    )
    class(erg) <- ("simex")
    if (class(model)[1] == "polr")
      erg$zeta <- z[1,-1][(length(coef(model))+1):length(z[1,-1])]
    type <- 'response'
    if (class(model)[1] == "polr")
      type <- 'probs'
    fitted.values <- predict(erg, newdata = model$model[, -1, drop = FALSE],
                             type = type)
    erg$fitted.values <- fitted.values
    if (class(model)[1] == "polr") {
      erg$residuals <- NULL
    } else if (is.factor(model$model[, 1])) {
      erg$residuals <-
        as.numeric(levels(model$model[, 1]))[model$model[, 1]] - fitted.values
    } else {
      erg$residuals <- model$model[, 1] - fitted.values
    }
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

#' @describeIn simex Plot the simulation and extrapolation step
#' @export
plot.simex <- function(x, xlab = expression((1 + lambda)), ylab = colnames(b[,
                                                                             -1]), ask = FALSE, show = rep(TRUE, NCOL(b) - 1), ...) {
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
         nonl = for (i in 1:length(p.names)) d[, i] <- predict(x$extrapolation[[p.names[i]]], newdata = data.frame(lambda = a)),
         log2 = for (i in 1:length(p.names)) d[, i] <- predict(x$extrapolation[[p.names[i]]], newdata = data.frame(lambda = a)) -
              ((abs(apply(x$SIMEX.estimates[-1, -1], 2, min)) + 1) * (apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0))[i],
         logl = d <- t(t(exp(predict(x$extrapolation, newdata = data.frame(lambda = a)))) - ((abs(apply(x$SIMEX.estimates[-1, -1], 2, min)) + 1) * (apply(x$SIMEX.estimates[-1, -1], 2, min) <= 0))))
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

#' @describeIn simex Predict using simex correction
#' @export
predict.simex <- function(object, newdata, ...) {
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


#' @describeIn simex Print simex nicely
#' @export
print.simex <- function(x, digits = max(3, getOption("digits") - 3), ...) {
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

#' @describeIn simex Print summary nicely
#' @export
print.summary.simex <- function(x, digits = max(3, getOption("digits") -
                                                  3), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNaive model: \n")
  print(x$naive.model)
  cat("\nSimex variable :\n")
  if (length(x$measurement.error) == length(x$SIMEXvariable)) {
    var.error <- t(matrix(x$measurement.error))
    dimnames(var.error) <- list("Measurement error :", x$SIMEXvariable)
    print(var.error)
  } else {
    m.e <- as.matrix(x$measurement.error)
    var.error <- matrix(ncol = NCOL(x$measurement.error))
    for (i in 1:NCOL(x$measurement.error)) {
      ifelse(length(unique(m.e[, i])) == 1, var.error[, i] <- unique(m.e[, i]), var.error[, i] <- "heteroscedastic")
    }
    dimnames(var.error) <- list("Measurement error :", x$SIMEXvariable)
    print(var.error)
  }
  cat("\n\nNumber of iterations: ", x$B, "\n")
  if (!is.null(x$residuals)) {
    cat("\nResiduals: \n")
    print(summary(x$residuals), digits)
  }
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


#' @describeIn simex Refits the model with a different extrapolation function
#' @export
refit.simex <- function(object, fitting.method = "quadratic", jackknife.estimation = "quadratic",
                        asymptotic = TRUE, ...)
        .refit(object, fitting.method = fitting.method,
             jackknife.estimation = jackknife.estimation, asymptotic = asymptotic,
             allowed.fitting = c("quad", "line", "nonl"), allowed.jackknife = c("quad", "line", "nonl", FALSE), ...)

#' @describeIn simex Summary of simulation and extrapolation
#' @export
summary.simex <- function(object, ...) {
  if (class(object$model)[1] == "polr")
    est <- c(coef(object), object$zeta)
  else
    est <- coef(object)
  p.names <- names(est)
  est.table <- list()
  if (is.null(resid(object)))
    n <- object$model$n
  else
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
  class(ans) <- "summary.simex"
  ans$coefficients <- est.table
  ans$residuals <- resid(object)
  ans$call <- object$call
  ans$B <- object$B
  ans$naive.model <- object$model$call
  ans$SIMEXvariable <- object$SIMEXvariable
  ans$measurement.error <- object$measurement.error
  return(ans)
}






