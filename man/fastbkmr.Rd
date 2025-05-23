% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastbkmr.R
\name{fastbkmr}
\alias{fastbkmr}
\title{Approximate Bayesian Kernel Machine Regression via Random Fourier Features}
\usage{
fastbkmr(
  y,
  X,
  Z,
  Z.eq = as.formula("~ .-1"),
  K = 20000,
  J = 20,
  a = 0.001,
  b = 0.001,
  s2_gamma = 100,
  theta = rep(0.5, ncol(X)),
  sig2 = 1,
  tau2 = 1,
  verbose = T,
  omega_tune_rate = 0.2,
  beta_tune_rate = 0.2
)
}
\arguments{
\item{y}{A numeric vector for the outcome.}

\item{X}{A numeric matrix or data.frame for exposure variables.}

\item{Z}{A numeric matrix or data.frame for covariates (confounders).}

\item{Z.eq}{A formula specifying covariate transformation (default: ~ . -1).}

\item{K}{Number of MCMC iterations (default: 20000).}

\item{J}{Number of basis functions (default: 20).}

\item{a}{Shape parameter for inverse-gamma priors (default: 0.001).}

\item{b}{Rate parameter for inverse-gamma priors (default: 0.001).}

\item{s2_gamma}{Prior variance for fixed effects (default: 100).}

\item{theta}{Initial theta values (default: 0.5).}

\item{sig2}{Initial sigma^2 (default: 1).}

\item{tau2}{nitial tau^2 (default: 1).}

\item{verbose}{Logical; whether to print tuning info (default: TRUE).}

\item{omega_tune_rate}{Tuning rate for omega (default: 0.2).}

\item{beta_tune_rate}{Tuning rate for beta (default: 0.2).}
}
\value{
A list with posterior summaries and samples.
}
\description{
Approximate Bayesian Kernel Machine Regression via Random Fourier Features
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
}
