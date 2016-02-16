## Auxiliary functions for use with signal::butter, the Butterworth filter

#' Initialize parameters of the pm_opt_filter on evaluation scale.
#' In principle, we could/should do something data based here?  Instead, I just made up
#' default numbers.
#' 
#' @param prev_phi list of previous filtering parameters for the given random variable
#' @param x vector of observations for the given random variable to be filtered
#' @param ... absorb arguments
#' 
#' @return initialized list of parameters
initialize_filter_params_butterworth_filter <- function(prev_phi,
    x,
    ...) {
    
    new_phi <- prev_phi
    
    default_w <- 0.3
    new_phi$logit_w <- log(default_w) - log(1 - default_w)
    new_phi$W <- default_w
    
    return(new_phi)
}

#' Vectorize the parameters of the pm_opt_filter and convert
#' to estimation scale.
#' 
#' @param phi list of parameters on evaluation scale
#' 
#' @return vector of parameters on estimation scale
vectorize_filter_params_butterworth_filter <- function(phi) {
    return(phi$logit_w)
}

#' Update the parameter list phi from a vector of parameters on estimation scale
#' 
#' @param phi_est_vector a numeric vector of length at least 2
#' @param phi a list of parameters to update
#' 
#' @return updated list phi
update_filter_params_from_vectorized_butterworth_filter <- function(phi_est_vector, phi) {
    ## update parameters on logit scale -- so that they can be returned
    ## easily by vectorize_filter_params_butterworth_filter
    phi$logit_w <- phi_est_vector[1]
    
    ## update on frequency scale: w = exp(logit_w) / [1 + exp(logit_w)]
    log_denom <- logspace_sum(c(0, phi$logit_w))
    phi$W <- exp(phi$logit_w - log_denom)
    
    return(list(
            phi = phi,
            num_phi_vals_used = 1L
        ))
}

#' Compute filter arguments for the Butterworth filter based on parameters phi
#' 
#' @param phi list of parameters including n, W
#' @param x vector to apply filtration to
#' 
#' @return named list with parameters for a call to signal::filter.  Filter
#'   coefficients are computed via signal::remez
compute_filter_args_butterworth_filter <- function(phi, x) {
    
    ## Obtain Butterworth filter coefficients
    ## If W is too small or too large, we can get numerical instability in the filtering.
    ## Check for "extreme" values and manually set filter coefficients in those cases.
    ## For W too extreme, we do "no filtering".  The cutoffs for W were determined experimentally.
    ## 0.05 corresponds roughly to the value at which MA coefficients go below machine epsilon.
    ## 0.95 corresponds roughly to the value at which numerical instability shows up in the
    ##  national flu data set.
    if(phi$W < 0.05 || phi$W > 0.95) {
        filt <- c(1, rep(0, phi$n - 1))
        class(filt) <- "Ma"
    } else {
        filt <- signal::butter(
            n = phi$n,
            W = phi$W
        )
    }
    
    return(list(
        x = x,
        filt = filt,
        impute_fn = phi$impute_fn
    ))
}
