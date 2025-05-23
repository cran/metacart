% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/re_find_cutoff.R
\name{re.cutoff_cpp}
\alias{re.cutoff_cpp}
\title{A function to find the split point}
\usage{
re.cutoff_cpp(g, vi, x, inx.s, cnode, minbucket)
}
\arguments{
\item{g}{the effect size}

\item{vi}{the sampling variance}

\item{x}{the splitting moderator}

\item{inx.s}{indicates whether a study belongs to the candidate parent leaf}

\item{cnode}{the terminal nodes that the studies belong to in the current tree}

\item{minbucket}{the minimum number of the studies in a terminal node}
}
\value{
a vector including the split point, Q, and tau2
}
\description{
A function to find the split point
}
\keyword{internal}
