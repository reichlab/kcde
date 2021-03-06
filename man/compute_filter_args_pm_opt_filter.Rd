% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parks-mcclellan-opt-filter.R
\name{compute_filter_args_pm_opt_filter}
\alias{compute_filter_args_pm_opt_filter}
\title{Compute filter arguments for the Parks-McClellan optimal FIR filter based on
parameters phi}
\usage{
compute_filter_args_pm_opt_filter(phi, x)
}
\arguments{
\item{phi}{list of parameters including n, f_c_0, and f_c_1}

\item{x}{vector to apply filtration to}
}
\value{
named list with parameters for a call to signal::filter.  Filter
  coefficients are computed via signal::remez
}
\description{
Compute filter arguments for the Parks-McClellan optimal FIR filter based on
parameters phi
}

