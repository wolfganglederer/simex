
extract.covmat <- function(model) {
  type <- class(model)[1]
  sum.model <- summary(model)
  switch(type , glm = covmat <- sum.model$cov.scaled
              , lm  = covmat <- sum.model$cov.unscaled * sum.model$sigma^2
              , gam = covmat <- model$Vp
              , nls = covmat <- sum.model$cov.unscaled * sum.model$sigma^2
              , polr= covmat <- stats::vcov(model)
              , lme = covmat <- model$apVar
              , nlme = covmat <- model$apVar)
  return(covmat)
}


fit.logl <- function(lambda, p.names, estimates) {
  est <- (t(t(estimates) + (abs(apply(estimates, 2, min)) + 1) * (apply(estimates,
                                                                        2, min) <= 0)))
  # calculating starting values
  start.mod <- stats::lm(I(log(t(t(estimates) + (abs(apply(estimates, 2, min)) + 1) * (apply(estimates, 2, min) <= 0)))) ~ lambda)
  start.val <- stats::coef(start.mod)
  extrapolation <- list()
  # doing the extrapolation step
  for (d in p.names) {
    extrapolation[[d]] <- try(stats::nls(est[, d] ~ exp(gamma.0 + (gamma.1 * lambda)), start = list(gamma.0 = start.val[1, d], gamma.1 = start.val[2, d])), silent = TRUE)
  }
  # security, in the case that nls() does not converge the 'logl'-method
  # is used
  if (any(lapply(extrapolation, class) == "try-error")) {
    warning("Function nls() did not converge, using 'logl' as fitting method",
            call. = FALSE)
    extrapolation <- start.mod
  }
  return(extrapolation)
}


fit.nls <- function (lambda, p.names, estimates){
  extrapolation <- list()
  # lambda values for the interpolation
  lambdastar <- c(0, max(lambda) / 2, max(lambda))
  for (d in p.names) {
    # calculating starting values
    # quadratic interpolation
    extrapolation.quad <- stats::lm(estimates[, d] ~ lambda + I(lambda^2))
    # interpolation for the values of lambdastar
    a.nls <- stats::predict(extrapolation.quad,
                     newdata = data.frame(lambda = lambdastar))
    # analytic 3-point fit => good starting values
    gamma.est.3 <- ((a.nls[2] - a.nls[3]) * lambdastar[3] *
                      (lambdastar[2] - lambdastar[1]) - lambdastar[1] *
                      (a.nls[1] - a.nls[2]) * (lambdastar[3] - lambdastar[2])) /
      ((a.nls[1] - a.nls[2]) * (lambdastar[3] - lambdastar[2]) -
         (a.nls[2] - a.nls[3]) * (lambdastar[2] - lambdastar[1]))
    gamma.est.2 <- ((a.nls[2] - a.nls[3]) * (gamma.est.3 + lambdastar[2]) *
                      (gamma.est.3 + lambdastar[3])) / (lambdastar[3] - lambdastar[2])
    gamma.est.1 <- a.nls[1] - (gamma.est.2 / (gamma.est.3 + lambdastar[1]))
    # fitting the nls-model for the various coefficients
    extrapolation[[d]] <-
      stats::nls(estimates[, d] ~ gamma.1 + gamma.2 / (gamma.3 + lambda),
          start = list(gamma.1 = gamma.est.1, gamma.2 = gamma.est.2,
                       gamma.3 = gamma.est.3))
  }
  return(extrapolation)
}


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
                                                 prob = mc[, i], replace = TRUE)
    }
  }
  return(data.mc)
}


construct.s <- function(ncoef, lambda, fitting.method, extrapolation = NULL) {
  nl <- length(lambda)
  switch(fitting.method, quad = ngamma <- 3, line = ngamma <- 2, nonl = ngamma <- 3,
         logl = ngamma <- 2, log2 = ngamma <- 2)
  # construct a matrix of 0
  null.mat <- matrix(0, nrow = ngamma, ncol = ncoef)
  s <- list()
  for (j in 1:ncoef) {
    # define the first matrix
    switch(fitting.method,
           quad = gamma.vec <- c(-1, -lambda[1], -lambda[1]^2),
           line = gamma.vec <- c(-1, -lambda[1]),
           nonl = gamma.vec <- c(1, 1/(coef(extrapolation[[j]])[3] + lambda[1]), -coef(extrapolation[[j]])[2]/(coef(extrapolation[[j]])[3] + lambda[1])^2),
           logl = gamma.vec <- c(exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] * lambda[1]), exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] * lambda[1]) * lambda[1]),
           log2 = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] + coef(extrapolation[[j]])[2] * lambda[1]), exp(coef(extrapolation[[j]])[1] + coef(extrapolation[[j]])[2] * lambda[1]) * lambda[1]))
    a <- null.mat
    # exchange the column j with the deviation of G(gamma, lambda)
    a[, j] <- gamma.vec
    for (i in 2:nl) {
      switch(fitting.method,
             quad = gamma.vec <- c(-1, -lambda[i], -lambda[i]^2),
             line = gamma.vec <- c(-1, -lambda[i]),
             nonl = gamma.vec <- c(1, 1/(coef(extrapolation[[j]])[3] + lambda[i]), -coef(extrapolation[[j]])[2]/(coef(extrapolation[[j]])[3] + lambda[i])^2),
             logl = gamma.vec <- c(exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] * lambda[i]), exp(coef(extrapolation)[1, j] + coef(extrapolation)[2, j] * lambda[i]) * lambda[i]),
             log2 = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] + coef(extrapolation[[j]])[2] * lambda[i]), exp(coef(extrapolation[[j]])[1] + coef(extrapolation[[j]])[2] * lambda[i]) * lambda[i]))
      b <- null.mat
      b[, j] <- gamma.vec
      a <- cbind(a, b)
    }
    s[[j]] <- a
  }
  s <- t(matrix(unlist(lapply(s, t), recursive = FALSE), nrow = nl * ncoef, ncol = ngamma * ncoef))
  return(s)
}



