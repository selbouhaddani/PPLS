% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossval_PPLS.R
\name{cv_PPLS}
\alias{cv_PPLS}
\title{Predicts X or Y with PPLS}
\usage{
cv_PPLS(X, Y, nr_comp, nr_folds, ...)
}
\arguments{
\item{X}{X data}

\item{Y}{Y data}

\item{nr_comp}{Strictly positive integer. The number of PPLS components to fit.}

\item{nr_folds}{Positive integer. Number of folds to consider. Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two folds does not make sense.}

\item{...}{PPLS arguments}
}
\value{
Cross-validated RMSEP averaged over all elements
}
\description{
Predicts X or Y based on new data on Y or X
}
