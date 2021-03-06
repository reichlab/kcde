% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{interior_linear_interpolation}
\alias{interior_linear_interpolation}
\title{Fill interior NAs in x by linear interpolation}
\usage{
interior_linear_interpolation(x)
}
\arguments{
\item{x}{a numeric vector of values to fill by linear interpolation}
}
\value{
a matrix with nrow = length(x) and one column, with interior NAs in x
  replaced by linearly interpolated values
}
\description{
Fill interior NAs in x by linear interpolation
}

