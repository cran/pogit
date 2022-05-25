\name{pogit-package}
\alias{pogit-package}
\docType{package}
\title{Bayesian variable selection for a Poisson-Logistic model}

\description{
This package provides Bayesian variable selection 
for regression models of under-reported count data as well as for (overdispersed) 
Poisson, negative binomial and binomial logit regression models using spike and slab priors. 
For posterior inference, MCMC sampling schemes are used that rely on data 
augmentation and/or auxiliary mixture sampling techniques. 
Details can be found in Dvorzak and Wagner (2016). 
%A concise (1-5 lines) description of the package 
}
\details{
%\tabular{ll}{
%Package: \tab pogit\cr
%Type: \tab Package\cr
%Version: \tab 1.0.0\cr
%Date: \tab 2015-04-17\cr
%License: \tab GPL-2\cr
%}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
The main function is \code{\link{pogitBvs}} which provides Bayesian variable
selection for a Poisson-Logistic (Pogit) model to account for potential 
under-reporting of count data. The Pogit model, introduced by Winkelmann 
and Zimmermann (1993), is specified by combining a Poisson model for the data 
generating process of counts and a logit model for the fallible reporting 
process, where the outcomes of both processes may depend on a set of 
potential covariates. 
By augmenting the observed data with the unobserved counts, the model 
can be factorized into a Poisson and a binomial logit model part. Hence, 
the MCMC sampling algorithm for this two-part model is based on  
data augmentation and sampling schemes for a Poisson and a binomial
logit model. 

Though part of the main function, the functions \code{\link{poissonBvs}} 
and \code{\link{logitBvs}} can be used separately to perform 
Bayesian variable selection for Poisson or binomial logit regression models. 
An alternative to \code{\link{poissonBvs}} is provided by the function 
\code{\link{negbinBvs}} to deal with overdispersion of count data.  
The sampling algorithms are based on auxiliary mixture sampling 
techniques. 


All functions return an object of class "\code{pogit}" with methods
\code{\link{print.pogit}}, \code{\link{summary.pogit}} and 
\code{\link{plot.pogit}} to summarize and display the results.  
}
\author{
Michaela Dvorzak <m.dvorzak@gmx.at>, Helga Wagner

Maintainer: Michaela Dvorzak <m.dvorzak@gmx.at>
}
\references{
Dvorzak, M. and Wagner, H. (2016). Sparse Bayesian modelling of underreported
count data. \emph{Statistical Modelling}, \strong{16}(1), 24 - 46, 
\doi{10.1177/1471082x15588398}.

Winkelmann, R. and Zimmermann, K. F. (1993). Poisson-Logistic regression. 
Department of Economics, University of Munich, Working Paper No. 93 - 18.
}

\keyword{ package }
\seealso{
\code{\link{pogitBvs}}, \code{\link{logitBvs}}, \code{\link{poissonBvs}}, 
\code{\link{negbinBvs}}
}
\examples{
## see examples for pogitBvs, logitBvs, poissonBvs and negbinBvs
}
