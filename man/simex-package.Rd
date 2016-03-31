\name{simex-package}
\Rdversion{1.1}
\encoding{latin1}
\alias{simex-package}
\docType{package}

\title{Error or misclassification correction in models using (MC)SIMEX}
\description{
  Package \code{simex} is an implementation of the SIMEX--algorithm by Cook
  and Stephanski and the MCSIMEX--Algorithm by Küchenhof, Mwalili and Lesaffre.

}

\details{
  \tabular{ll}{
    Package: \tab simex\cr
    Type: \tab Package\cr
    Version: \tab 1.7\cr
    Date: \tab 2016-03-25\cr
    License: \tab GPL 2 or above\cr
    LazyLoad: \tab yes\cr
  }
  The package includes first of all the implementation for the SIMEX-- and
  MCSIMEX--Algorithms. Jackknife and asymptotic variance estimation are
  implemented. Various methods and analytic tools are provided for a simple and
  fast access to the SIMEX-- and MCSIMEX--Algorithm.

  Functions \code{simex()} and \code{mcsimex()} can be used on models issued
  from \code{lm()}, \code{glm()} with asymtotic estimation.
  Models from \code{nls()}, \code{gam()} (package \pkg{mgcv}),
  \code{polr()} (package \pkg{MASS}),
  \code{lme()} and \code{nlme()} (package \pkg{nlme}) can also be corrected with
  these algorithms, but without asymptotic estimations.
}

\author{
Wolfgang Lederer, Heidi Seibold, Helmut Küchenhoff

Maintainer: Wolfgang Lederer,\email{wolfgang.lederer@gmail.com}
}

\references{
Lederer, W. and Küchenhoof, H. (2006) A short introduction to the SIMEX and
  MCSIMEX. \emph{R News}, \bold{6/4}, 26 -- 31
}

\seealso{
  \code{\link[simex]{simex}}, \code{\link[simex]{mcsimex}}, \code{\link[simex]{misclass}}

  and for functions generating the initial naive models:
  \code{\link[stats]{lm}}, \code{\link[stats]{glm}}, \code{nls},
  \code{gam}, \code{lme}, \code{nlme}, \code{polr}
}

\examples{
# See example(simex) and example(mcsimex)
## Seed
set.seed(49494)

## simulating the measurement error standard deviations
sd_me1 <- 0.3
sd_me2 <- 0.4
temp <- runif(100, min = 0, max = 0.6)
sd_me_het1 <- sort(temp)
temp2 <- rnorm(100, sd = 0.1)
sd_me_het2 <- abs(sd_me_het1 + temp2)

## simulating the independent variables x (real and with measurement error):
x_real1 <- rnorm(100)
x_real2 <- rpois(100, lambda = 2)
x_real3 <- -4*x_real1 + runif(100, min = -2, max = 2)  # correlated to x_real

x_measured1 <- x_real1 + sd_me1 * rnorm(100)
x_measured2 <- x_real2 + sd_me2 * rnorm(100)
x_het1 <- x_real1 + sd_me_het1 * rnorm(100)
x_het2 <- x_real3 + sd_me_het2 * rnorm(100)

## calculating dependent variable y:
y1  <- x_real1 + rnorm(100, sd = 0.05)
y2 <- x_real1 + 2*x_real2 + rnorm(100, sd = 0.08)
y3 <- x_real1 + 2*x_real3 + rnorm(100, sd = 0.08)


### one variable with homoscedastic measurement error
(model_real <- lm(y1  ~ x_real1))

(model_naiv <- lm(y1  ~ x_measured1, x = TRUE))

(model_simex <- simex(model_naiv, SIMEXvariable = "x_measured1", measurement.error = sd_me1))
plot(model_simex)


### two variables with homoscedastic measurement errors
(model_real2 <- lm(y2 ~ x_real1 + x_real2))

(model_naiv2 <- lm(y2 ~ x_measured1 + x_measured2, x = TRUE))

(model_simex2 <- simex(model_naiv2, SIMEXvariable = c("x_measured1", "x_measured2"), 
                       measurement.error = cbind(sd_me1, sd_me2)))

plot(model_simex2)


### one variable with increasing heteroscedastic measurement error
model_real 

(mod_naiv1 <- lm(y1  ~ x_het1, x = TRUE))

(mod_simex1 <- simex(mod_naiv1, SIMEXvariable = "x_het1",
     measurement.error = sd_me_het1, asymptotic = FALSE))

plot(mod_simex1)

\dontrun{
### two correlated variables with heteroscedastic measurement errors
(model_real3 <- lm(y3 ~ x_real1 + x_real3))

(mod_naiv2 <- lm(y3 ~ x_het1 + x_het2, x = TRUE))

(mod_simex2 <- simex(mod_naiv2, SIMEXvariable = c("x_het1", "x_het2"), 
                     measurement.error = cbind(sd_me_het1, sd_me_het2), asymptotic = FALSE))
plot(mod_simex2)


### two variables, one with homoscedastic, one with heteroscedastic measurement error
model_real2

(mod_naiv3 <- lm(y2 ~ x_measured1 + x_het2, x = TRUE))

(mod_simex3 <- simex(mod_naiv3, SIMEXvariable = c("x_measured1", "x_het2"), 
                     measurement.error = cbind(sd_me1, sd_me_het2), asymptotic = FALSE))


### glm: two variables, one with homoscedastic, one with heteroscedastic measurement error
t <- x_real1 + 2*x_real2 
g <- 1 / (1 + exp(-t))
u <- runif(100)
ybin <- as.numeric(u < g)


(logit_real <- glm(ybin ~ x_real1 + x_real2, family = binomial))

(logit_naiv <- glm(ybin ~ x_measured1 + x_het2, x = TRUE, family = binomial))

(logit_simex <- simex(logit_naiv, SIMEXvariable = c("x_measured1", "x_het2"), 
                      measurement.error = cbind(sd_me1, sd_me_het2), asymptotic = FALSE))
summary(logit_simex)
print(logit_simex)
plot(logit_simex)
}
}

\keyword{package}
\keyword{models}
\keyword{regression}
