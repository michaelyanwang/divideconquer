% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestep.R
\name{iteration.fun}
\alias{iteration.fun}
\title{One-step iteration function for divide-and-conquer logistic model}
\usage{
iteration.fun(dat.list, bini, kk.list)
}
\arguments{
\item{dat.list}{list of subsets (after dividing). In each subset, first column = outcome, rest = design matrix}

\item{bini}{initial estimator as starting point}

\item{kk.list}{which subsets of dat.list get one-step update}
}
\value{
a list with b.k a matrix of one step estimator and Ahat the negative information matrix
}
\description{
\code{iteration.fun} estimates a one-step for a logistic regression for each subset.
}
\examples{
iteration.fun(dat.list=dat.list,bini=bini,kk.list=2:K)
}
\author{
Yan Wang, Tianxi Cai
}
