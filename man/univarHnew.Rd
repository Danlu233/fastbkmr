% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univarHnew.R
\name{univarHnew}
\alias{univarHnew}
\title{Estimate univariate exposure-response function on new grid of points}
\usage{
univarHnew(
  fit,
  whichX,
  Xgrid = NULL,
  fixlevel = 0.5,
  center = NULL,
  ngrid = 50,
  alpha = 0.05
)
}
\arguments{
\item{fit}{An object containing the results from \code{fastbkmr} function.}

\item{whichX}{Vector identifying which exposure should be selected to varying. The number of selected exposures should be one.}

\item{Xgrid}{Values at grid points to cover the range of selected exposure.}

\item{fixlevel}{Quantile that other exposures will be fixed at (default: 0.5).}

\item{center}{The quantile of the selected exposure at which the exposure–response function equals zero.}

\item{ngrid}{The number of grids points to cover the range of selected exposure (default: 50).}

\item{alpha}{100(1-alpha)\% posterior interval (default: 0.05).}
}
\value{
A data frame with selected exposure name, exposure value, posterior mean estimate, and 100(1-alpha)\% posterior interval.
}
\description{
Estimate univariate exposure-response function on new grid of points
}
\examples{
data(simudat)
fit <- fastbkmr(y = simudat$y,
X = simudat[,1:5],
Z = simudat[,6:10],
Z.eq = as.formula("~ z1 + as.factor(z2) + z3 + z4 + as.factor(z5)"),
K = 1000,
J = 10,
verbose = FALSE)
est_uni <- univarHnew(fit = fit, whichX = "x1")
}
