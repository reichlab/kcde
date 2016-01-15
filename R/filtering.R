## Functions for interaction between kcde estimation and prediction routines and linear filtering
## 
## initialize_phi
## extract_vectorized_phi_est_from_phi
## update_phi_from_vectorized_phi_est
## apply_filters

#' Initialize parameter values for filtering
#'
#' @param prev_phi previous phi values before update to variables and lags
#' @param updated_vars_and_offsets lags after update
#' @param update_var_name character; the name of the variable added/removed from the model
#' @param update_offset_value integer; the offset that was added/removed from the model
#' @param update_offset_type integer; the type of the offset that was added/removed from the model
#' @param data data matrix
#' @param kcde_control list of control parameters for kcde
#' 
#' @return list of phi parameters
initialize_phi <- function(prev_phi,
    updated_vars_and_offsets,
    update_var_name,
    update_offset_value,
    update_offset_type,
    data,
    kcde_control) {
    
    phi <- lapply(seq_along(kcde_control$filter_control),
        function(ind) {
            updated_var_in_current_component <-
                any(update_var_name == paste0("filtered_", kcde_control$filter_control[[ind]]$var_name) &
                        update_offset_type == "lag")
            if(updated_var_in_current_component) {
                potential_cols_in_component <-
                    kcde_control$filter_control[[ind]]$var_name
                cols_used <- colnames(data) %in% potential_cols_in_component
                if(any(cols_used)) {
                    fn_name <- kcde_control$filter_control[[ind]]$
                        initialize_filter_params_fn
                    
                    fn_args <- kcde_control$filter_control[[ind]]$
                        initialize_filter_params_args
                    fn_args$update_var_name <- update_var_name
                    fn_args$update_offset_value <- update_offset_value
                    fn_args$update_offset_type <- update_offset_type
                    fn_args$kcde_control <- kcde_control
                    fn_args$x <- data[, cols_used]
                    
                    ## the business with temp_prev_phi here is necessary
                    ## because on the first call for a given filter component,
                    ## prev_phi may be NULL, but on later calls, prev_phi
                    ## will have already been initialized with the values in
                    ## phi_fixed
                    temp_prev_phi <- prev_phi[[ind]]
                    temp_prev_phi[
                        names(kcde_control$filter_control[[ind]]$fixed_filter_params)
                    ] <- kcde_control$filter_control[[ind]]$fixed_filter_params
                    fn_args$prev_phi <- temp_prev_phi
                    
                    return(do.call(fn_name, fn_args))
                } else {
                    return(NULL)
                }
            } else {
                return(prev_phi[[ind]])
            }
        }
    )
    
    return(phi)
}


#' Extract a vector of parameter values that are to be estimated from phi,
#' on the estimation scale.
#' 
#' @param phi filtering parameters phi in list form
#' @param filter_control control parameters for the kcde fit
#' 
#' @return numeric vector with parameter values
extract_vectorized_phi_est_from_phi <- function(phi,
    filter_control) {
    
    phi_vector <- c()
    
    for(ind in seq_along(filter_control)) {
        ## parameters that are being estimated
        if(!is.null(phi[[ind]])) {
            phi_vector <- c(phi_vector, do.call(
                    filter_control[[ind]]$vectorize_filter_params_fn,
                    c(list(phi = phi[[ind]]),
                        filter_control[[ind]]$vectorize_filter_params_args)
                ))
        }
    }
    
    return(phi_vector)
}

#' Convert phi from vector form to list form
#' 
#' @param phi_est_vector numeric vector of filter parameters phi that are
#'     being estimated, on estimation scale.
#' @param kcde_control control parameters for the kcde fit
#' 
#' @return list of lists of parameter values -- outer list has one component
#'     for each filter function, inner list has one component
#'     for each parameter used in the corresponding filter function
update_phi_from_vectorized_phi_est <- function(phi_est_vector,
	phi,
    filter_control) {
    
    phi_vector_component_start_ind <- 1L
    for(ind in seq_along(filter_control)) {
        ## parameters that are being estimated
        if(!is.null(phi[[ind]])) {
            temp <- do.call(
                filter_control[[ind]]$
                    update_filter_params_from_vectorized_fn,
                c(list(phi_est_vector = phi_est_vector[
                            seq(from = phi_vector_component_start_ind,
                                to = length(phi_est_vector))],
                        phi = phi[[ind]]),
                    filter_control[[ind]]$
                        update_filter_params_from_vectorized_args)
            )
            
            phi_vector_component_start_ind <- phi_vector_component_start_ind +
                temp$num_phi_vals_used
            
            phi[[ind]] <- temp$phi
        }
    }

    return(list(phi = phi,
        next_param_ind = phi_vector_component_start_ind))
}


#' Compute filtered values
#' 
#' @param data a data frame with original data (no filtering done)
#' @param phi a list with one component for each component of the filter
#'     function.  This component is a named list with arguments to the
#'     corresponding filter function.
#' @param kcde_control a list of kcde_control parameters for kcde
#' @param log boolean; if TRUE (default), return filter values on the log scale
compute_filter_values <- function(data,
    filter_control,
    phi) {
    
    ## Determine for which variables we need to perform filtering
    non_null_phi_inds <- which(!sapply(phi, is.null))
    var_names_to_filter <- sapply(non_null_phi_inds, function(ind) {
        filter_control[[ind]]$var_name
    })
    
    ## Allocate matrix to store filtered variables
    filtered_data <- matrix(NA, nrow = nrow(data), ncol = length(var_names_to_filter))
    colnames(filtered_data) <- paste0("filtered_", var_names_to_filter)
    
    ## Fill in filtered_data
    for(ind in seq_along(non_null_phi_inds)) {
        ## assemble arguments to filter coef function
        filter_args_args <- list(phi = phi[[non_null_phi_inds[ind]]])
        filter_args_args$x <- data[, var_names_to_filter[ind]]
        
        ## call filter coef function
        filter_args <-
            do.call(filter_control[[non_null_phi_inds[ind]]]$filter_args_fn,
                filter_args_args)
        
        ## do filtering and store result in filtered_data
        filtered_data[, non_null_phi_inds[ind]] <-
            do.call(filter_control[[non_null_phi_inds[ind]]]$filter_fn, filter_args)
    }
    
    ## return
    return(cbind(data, as.data.frame(filtered_data)))
}

#' Forward and reverse pass filter to correct phase shift introduced by one-pass filter.
#' Closely based on signal::filtfilt -- see the documentation there for further discussion.
#' Two main differences with the implementation in signal::filtfilt:
#'   (1) We pad the end of the time series with the last observed value instead
#'       of zeros.  Based on very informal plotting, this seems to give better
#'       values at the end of the time series, which will be important for
#'       prediction
#'   (2) We call stats::filter to do the filtering instead of signal::filter.
#'       This is because signal::filter has an enforced call to na.omit, which
#'       clashes with our handling of NA values.  On the other hand, I think
#'       we should do something like linear or quadratic interpolation of
#'       internal NAs.  This would solve the problem with calling
#'       signal::filter, which is preferred since it is more flexible.
#' 
#' @param filter vector of filter coefficients for a FIR filter
#' @param x vector of data to be filtered
#' @param method method for a call to stats::filter -- see the documentation at stats::filter
#' @param sides sides for a call to stats::filter -- see the documentation at stats::filter
#' @param circular circular for a call to stats::filter -- see the documentation at stats::filter
#' 
two_pass_filter <- function(filter, x, method, sides, circular) {
#    y <- stats::filter(filter = filter, x = c(x, numeric(2 * length(filter))), method = method, sides = sides, circular = circular)
    y <- stats::filter(filter = filter, x = c(x, rep(x[length(x)], 2 * length(filter))), method = method, sides = sides, circular = circular)
    y <- rev(stats::filter(filter = filter, x = rev(y), method = method, sides = sides, circular = circular))[seq_along(x)]
    return(y)
}
