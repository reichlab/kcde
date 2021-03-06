% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{kcde_crossval_estimate_parameter_loss}
\alias{kcde_crossval_estimate_parameter_loss}
\title{Using cross-validation, estimate the loss associated with a particular set
of lags and kernel function parameters.  We calculate log-score loss}
\usage{
kcde_crossval_estimate_parameter_loss(combined_params_vector, phi, theta,
  vars_and_offsets, data, leading_rows_to_drop, trailing_rows_to_drop,
  additional_rows_to_drop, kcde_control)
}
\arguments{
\item{combined_params_vector}{vector of parameters for filtering and kernel
functions that are being estimated}

\item{theta}{list of kernel function parameters, both those that are being
estimated and those that are out of date.  Possibly the values of
parameters being estimated are out of date; they will be replaced with
the values in theta_est_vector.}

\item{vars_and_offsets}{list representing combinations of variables and lags
included in the model}

\item{data}{the data frame to use in performing cross validation}

\item{kcde_control}{a list of parameters specifying how the fitting is done}
}
\value{
numeric -- cross-validation estimate of loss associated with the
    specified parameters
}
\description{
Using cross-validation, estimate the loss associated with a particular set
of lags and kernel function parameters.  We calculate log-score loss
}

