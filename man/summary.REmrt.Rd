% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.REmrt}
\alias{summary.REmrt}
\title{Summary of the results of a RE meta-tree object}
\usage{
\method{summary}{REmrt}(object, digits = 3, ...)
}
\arguments{
\item{object}{fitted tree of class \code{REmrt}.}

\item{digits}{specified number of decimals in the printed results.}

\item{\dots}{additional arguments to be passed.}
}
\value{
Summary of a RE meta-tree
}
\description{
Summary of the results of a RE meta-tree object
}
\details{
If no moderator effect is detected,
the summary function will show the standard meta-analysis results.
Otherwise, the summary function will show the subgroup meta-analysis results,
with the significance test results for moderator effects, the split points of the moderators,
and the estimated subgroup summary effect sizes.
}
