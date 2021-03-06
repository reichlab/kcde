% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel-interaction.R
\name{get_theta_optim_bounds}
\alias{get_theta_optim_bounds}
\title{Get vectors of lower and upper bounds for theta parameters to be used in
calls to optim.}
\usage{
get_theta_optim_bounds(theta, kcde_control)
}
\arguments{
\item{theta}{list of parameters theta that are being estimated}

\item{kcde_control}{list of control parameters for kcde}
}
\value{
list with two components: lower and upper, giving vectors
  with lower and upper bounds for possible parameter values
  Initialize parameter values.
}
\description{
Get vectors of lower and upper bounds for theta parameters to be used in
calls to optim.
}

