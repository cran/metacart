% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootstrapCIfunctions.R
\name{BootCI}
\alias{BootCI}
\title{This functions performs bootstrap to compute the confidence intervals for the subgroup effect size estimates.
This function is only applicable to Random effects metaregression trees with 2 terminal nodes or more.}
\usage{
BootCI(Metatree, nboot = 50)
}
\arguments{
\item{Metatree}{fitted tree of class \code{REmrt}.}

\item{nboot}{number of bootstrap samples.}
}
\value{
tree containing the input tree, the Bootstrap estimates for the effect sizes and standard errors, Bootstrap estimate for tau2, and the Bootstrap bias correction.
}
\description{
This functions performs bootstrap to compute the confidence intervals for the subgroup effect size estimates.
This function is only applicable to Random effects metaregression trees with 2 terminal nodes or more.
}
\examples{
set.seed(12345) 
data(dat.BCT2009)
library(Rcpp)
REtree <- REmrt(g ~ T1 + T2+ T4 +T25, vi = vi, data = dat.BCT2009, c.pruning = 0)
BootTree<-BootCI(REtree, nboot = 3)
summary(BootTree)

}
