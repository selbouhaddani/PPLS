% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_W_multi.R
\name{scores.PPLS}
\alias{scores.PPLS}
\title{Get scores from PPLS fit.}
\usage{
scores.PPLS(fit, X, Y, subset = NULL)
}
\arguments{
\item{fit}{A PPLS fit (an object of class PPLS)}

\item{X}{Matrix}

\item{Y}{Matrix}

\item{subset}{vector of positive integers denoting with components you want.}
}
\value{
Both the X and Y scores concatenated.
}
\description{
Get scores from PPLS fit.
}
\details{
If you want just the X or Y scores, subset the result i.e. \code{scores.PPLS(fit,X,Y)[1:nrow(X),]} or \code{scores.PPLS(fit,X,Y)[-(1:nrow(X)),]}.
}
