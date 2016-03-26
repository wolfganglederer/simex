# The simex Package

This pakage for R implements the simex procedure developed by Cook & Stefanski for dealing with measurement error models, as well as the mcsimex for misclassified data developed by Küchenhoff, Mwalili and Lesaffre.

It can be found on [CRAN] (https://cran.r-project.org/web/packages/simex/index.html), or the most current version can be installed via devtools: devtools::install_github("wolfganglederer/simex")


# simex NEWS:

## Version 1.7:
(by Wolfgang Lederer)

- Option fitting.method in function simex() now really works properly for "linear" (Wolfgang)
- Smaller documentation fixes (Wolfgang)
- Uploaded to Github and moved documentation to roxygen2 (Wolfgang)

## Version 1.5:
(by Heidi Seibold & Wolfgang Lederer)

- Option fitting.method in function simex() now works properly for "linear" (Wolfgang)
- Measurement errors may now be heteroscedastic. Therefore the input-type had to be changed. The measurement error now has to be a  matrix. (Heidi)
- The functions simex and print.summary.simex were adjusted. (Heidi)

## Version 1.4:
(by Ph. Grosjean <phgrosjean@sciviews.org>)

- Object classes were renamed 'simex' and 'mcsimex' to match their constructor's names simex() and mcsimex().
- Interface of mcsimex() has been homogenized with the one of simex(), with fitting.method at the sixth place instead of last one.
- print() methods now return x invisibly, as they are supposed to do.
- predict() methods failed when newdata was not provided. Fixed.
- refit() is now a generic function with methods for objects 'simex' and 'mcsimex'. Its arguments have been reworked to match arguments of simex() and mcsimex() functions (jackknife becomes jackknife.estimation).
- A NAMESPACE is added and functions that are not supposed to be used by the end-user are now hidden in the namespace (construct.s(), fit.log() and fit.nls()).
- Documentation has been rewritten to match current R standards. Methods are now documented in the same page as the object creator. 'Overview' is rewritten and renamed 'simex-package'.

# References

Küchenhoff, H., Mwalili, S. M.  and Lesaffre, E. (2006) A general method for dealing with misclassification in regression: The Misclassification SIMEX. *Biometrics*, **62**, 85 -- 96

Küchenhoff, H., Lederer, W. and E. Lesaffre. (2006) Asymptotic Variance Estimation for the Misclassification SIMEX. *Computational Statistics and Data Analysis*, **51**, 6197 -- 6211

Lederer, W. and Küchenhoff, H. (2006) A short introduction to the SIMEX and MCSIMEX. *R News*, **6(4)**, 26--31

Cook, J.R. and Stefanski, L.A. (1994) Simulation-extrapolation estimation in parametric measurement error models. *Journal of the American Statistical Association*, **89**, 1314 -- 1328

Carroll, R.J., Küchenhoff, H., Lombard, F. and Stefanski L.A. (1996) Asymptotics for the SIMEX estimator in nonlinear measurement error models. *Journal of the American Statistical Association*, **91**, 242 -- 250

Carrol, R.J., Ruppert, D., Stefanski, L.A. and Crainiceanu, C. (2006). *Measurement error in nonlinear models: A modern perspective.*, Second Edition. London: Chapman and Hall.
