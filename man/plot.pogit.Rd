% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.pogit}
\alias{plot.pogit}
\title{Plot an object of class \code{pogit}}
\usage{
\method{plot}{pogit}(
  x,
  type = "traceplot",
  burnin = TRUE,
  thin = TRUE,
  lag.max = NULL,
  ci = TRUE,
  maxPlots = NULL,
  ...
)
}
\arguments{
\item{x}{an object of class \code{pogit}}

\item{type}{type of plot: "\code{traceplot}" (default) for traceplots of the 
MCMC draws, "\code{acf}" for autocorrelation plots of the MCMC draws, 
"\code{density}" for density plots and "\code{hpd}" to display 
(model averaged) posterior means with 95\%-HPD intervals for the regression 
effects.}

\item{burnin}{logical. If \code{TRUE} (default), burn-in draws (as specified 
in \code{x}) are discarded.}

\item{thin}{logical. If \code{TRUE} (default), thinning (as specified in 
\code{x}) is considered for diagnostic MCMC plots.}

\item{lag.max}{maximum lag for autocorrelation plot; if \code{NULL} (default), 
the default of \code{\link[stats]{acf}} is used.}

\item{ci}{logical. If \code{TRUE} (default), the confidence interval in the
autocorrelation plot is shown (see \code{\link[stats]{acf}} for details).}

\item{maxPlots}{maximum number of plots on a single page; if \code{NULL} (default),
the number of plots dispayed on a single page is specified according to the
used model.}

\item{...}{further arguments (not used)}
}
\description{
This function provides traceplots, autocorrelation plots and density plots 
of the MCMC samples for an object of class "\code{pogit}" to graphically assess
convergence of the MCMC simulations. It also displays the (model averaged) 
posterior means and 95\%-HPD intervals for the regression effects.
}
\examples{
## see examples for pogitBvs, logitBvs, poissonBvs and negbinBvs
}
\author{
Michaela Dvorzak <m.dvorzak@gmx.at>
}
