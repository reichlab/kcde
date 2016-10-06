## Auxiliary functions for use with signal::remez, the Parks-McClellan optimal
## FIR filter.  This method is not useful for this package; the C implementation
## of the method in the signal package has some issues that we run into when
## attempting to use numerical methods to optimize band parameters.
##
## initialize_filter_params_pm_opt_filter
## vectorize_filter_params_pm_opt_filter
## update_filter_params_from_vectorized_pm_opt_filter
## compute_filter_args_pm_opt_filter

#' Initialize parameters of the pm_opt_filter on evaluation scale.
#' In principle, we could/should do something data based here?  Instead, I just made up
#' default numbers.
#'
#' @param prev_phi list of previous filtering parameters for the given random variable
#' @param x vector of observations for the given random variable to be filtered
#' @param ... absorb arguments
#'
#' @return initialized list of parameters
initialize_filter_params_pm_opt_filter <- function(prev_phi,
                                                   x,
                                                   ...) {
    new_phi <- prev_phi

    default_f_c_0 <- 0.2
    default_f_c_1 <- 0.3
    default_diff <- default_f_c_1 - default_f_c_0

    multilogit_denom <- log(1 - default_f_c_0 - default_diff)
    new_phi$multilogit_f_c_0 <-
        log(default_f_c_0) - multilogit_denom
    new_phi$multilogit_f_c_1_minus_f_c_0 = log(default_diff) - multilogit_denom

    new_phi$f_c_0 = default_f_c_0
    new_phi$f_c_1 = default_f_c_1

    return(new_phi)
}

#' Vectorize the parameters of the pm_opt_filter and convert
#' to estimation scale.
#'
#' @param phi list of parameters on evaluation scale
#'
#' @return vector of parameters on estimation scale
vectorize_filter_params_pm_opt_filter <- function(phi) {
    return(c(phi$multilogit_f_c_0, phi$multilogit_f_c_1_minus_f_c_0))
}

#' Update the parameter list phi from a vector of parameters on estimation scale
#'
#' @param phi_est_vector a numeric vector of length at least 2
#' @param phi a list of parameters to update
#'
#' @return updated list phi
update_filter_params_from_vectorized_pm_opt_filter <-
    function(phi_est_vector, phi) {
        ## update parameters on multilogit scale -- so that they can be returned
        ## easily by vectorize_filter_params_pm_opt_filter
        phi$multilogit_f_c_0 <- phi_est_vector[1]
        phi$multilogit_f_c_1_minus_f_c_0 <- phi_est_vector[2]

        ## update on frequency scale.  The relationships are:
        ## denom = [1 + exp(multilogit_f_c_0) + exp(multilogit_f_c_1_minus_f_c_0)]
        ## f_c_0 = exp(multilogit_f_c_0) / denom
        ## b = exp(multilogit_f_c_1_minus_f_c_0) / denom
        ## f_c_1 = f_c_0 + b
        log_denom <-
            logspace_sum(c(
                0,
                phi$multilogit_f_c_0,
                phi$multilogit_f_c_1_minus_f_c_0
            ))
        phi$f_c_0 <- exp(phi$multilogit_f_c_0 - log_denom)

        difference <- exp(phi$multilogit_f_c_1_minus_f_c_0 - log_denom)
        ## enforce minimum difference to avoid errors in signal::remez
        difference <- max(difference, 0.0001)
        phi$f_c_1 <- phi$f_c_0 + difference
        if (phi$f_c_1 > 1) {
            ## this may occur if we manually altered difference above
            phi$f_c_1 <- 1
            phi$f_c_0 <- phi$f_c_1 - difference
        }

        return(list(phi = phi,
                    num_phi_vals_used = 2L))
    }

#' Compute filter arguments for the Parks-McClellan optimal FIR filter based on
#' parameters phi
#'
#' @param phi list of parameters including n, f_c_0, and f_c_1
#' @param x vector to apply filtration to
#'
#' @return named list with parameters for a call to signal::filter.  Filter
#'   coefficients are computed via signal::remez
compute_filter_args_pm_opt_filter <- function(phi, x) {
    return(
        list(
            x = x,
            filter = signal::remez(
                n = phi$n,
                f = c(0, phi$f_c_0, phi$f_c_1, 1),
                a = c(1, 1, 0, 0),
                ftype = phi$ftype,
                density = phi$density
            ),
            method = "convolution",
            sides = 1L,
            circular = FALSE
        )
    )
}
