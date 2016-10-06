## Functions for interaction between kcde estimation and prediction routines and linear filtering
##
## initialize_phi
## extract_vectorized_phi_est_from_phi
## update_phi_from_vectorized_phi_est
## get_phi_optim_bounds
## compute_filter_values
## two_pass_stats_filter
## two_pass_signal_filter
## one_pass_signal_filter
## interior_linear_interpolation


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
                          any(
                              update_var_name == paste0("filtered_", kcde_control$filter_control[[ind]]$var_name) &
                                  update_offset_type == "lag"
                          )
                      if (updated_var_in_current_component) {
                          potential_cols_in_component <-
                              kcde_control$filter_control[[ind]]$var_name
                          cols_used <-
                              colnames(data) %in% potential_cols_in_component
                          if (any(cols_used)) {
                              fn_name <- kcde_control$filter_control[[ind]]$initialize_filter_params_fn

                              fn_args <- kcde_control$filter_control[[ind]]$initialize_filter_params_args
                              fn_args$update_var_name <- update_var_name
                              fn_args$update_offset_value <-
                                  update_offset_value
                              fn_args$update_offset_type <- update_offset_type
                              fn_args$kcde_control <- kcde_control
                              fn_args$x <- data[, cols_used]

                              ## the business with temp_prev_phi here is necessary
                              ## because on the first call for a given filter component,
                              ## prev_phi may be NULL, but on later calls, prev_phi
                              ## will have already been initialized with the values in
                              ## phi_fixed
                              temp_prev_phi <- prev_phi[[ind]]
                              temp_prev_phi[names(kcde_control$filter_control[[ind]]$fixed_filter_params)] <-
                                  kcde_control$filter_control[[ind]]$fixed_filter_params
                              fn_args$prev_phi <- temp_prev_phi

                              return(do.call(fn_name, fn_args))
                          } else {
                              return(NULL)
                          }
                      } else {
                          return(prev_phi[[ind]])
                      }
                  })

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

    for (ind in seq_along(filter_control)) {
        ## parameters that are being estimated
        if (!is.null(phi[[ind]])) {
            phi_vector <- c(phi_vector,
                            do.call(
                                filter_control[[ind]]$vectorize_filter_params_fn,
                                c(
                                    list(phi = phi[[ind]]),
                                    filter_control[[ind]]$vectorize_filter_params_args
                                )
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
    for (ind in seq_along(filter_control)) {
        ## parameters that are being estimated
        if (!is.null(phi[[ind]])) {
            temp <- do.call(
                filter_control[[ind]]$update_filter_params_from_vectorized_fn,
                c(
                    list(phi_est_vector = phi_est_vector[seq(from = phi_vector_component_start_ind,
                                                             to = length(phi_est_vector))],
                         phi = phi[[ind]]),
                    filter_control[[ind]]$update_filter_params_from_vectorized_args
                )
            )

            phi_vector_component_start_ind <-
                phi_vector_component_start_ind +
                temp$num_phi_vals_used

            phi[[ind]] <- temp$phi
        }
    }

    return(list(phi = phi,
                next_param_ind = phi_vector_component_start_ind))
}

#' Get vectors of lower and upper bounds for phi parameters to be used in
#' calls to optim.
#'
#' @param phi list of filter parameters phi that are being estimated
#' @param filter_control control parameters for filtering
#'
#' @return list with two components: lower and upper, giving vectors
#'   with lower and upper bounds for possible parameter values
get_phi_optim_bounds <- function(phi,
                                 filter_control) {
    lower <- NULL
    upper <- NULL

    for (ind in seq_along(filter_control)) {
        ## parameters that are being estimated
        if (!is.null(phi[[ind]])) {
            temp <- do.call(
                filter_control[[ind]]$get_phi_optim_bounds_fn,
                c(
                    list(phi = phi[[ind]]),
                    filter_control[[ind]]$get_phi_optim_bounds_args
                )
            )

            lower <- c(lower, temp$lower)
            upper <- c(upper, temp$upper)
        }
    }

    return(list(upper = upper,
                lower = lower))
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

    if (length(non_null_phi_inds) == 0) {
        ## No variables need filtering
        return(data)
    } else {
        ## At least one variable needs filtering

        var_names_to_filter <-
            sapply(non_null_phi_inds, function(ind) {
                filter_control[[ind]]$var_name
            })

        ## Allocate matrix to store filtered variables
        filtered_data <-
            matrix(NA,
                   nrow = nrow(data),
                   ncol = length(var_names_to_filter))
        colnames(filtered_data) <-
            paste0("filtered_", var_names_to_filter)

        ## Fill in filtered_data
        for (ind in seq_along(non_null_phi_inds)) {
            ## transform data to scale on which filtering will be performed
            if (!is.null(filter_control[[non_null_phi_inds[ind]]]$transform_fn)) {
                transformed_data_one_var <-
                    filter_control[[non_null_phi_inds[ind]]]$transform_fn(data[, var_names_to_filter[ind]])
            }

            ## assemble arguments to filter coef function
            filter_args_args <-
                list(phi = phi[[non_null_phi_inds[ind]]])
            filter_args_args$x <- transformed_data_one_var

            ## call filter coef function
            filter_args <-
                do.call(filter_control[[non_null_phi_inds[ind]]]$filter_args_fn,
                        filter_args_args)

            ## do filtering
            filtered_data_one_var <-
                do.call(filter_control[[non_null_phi_inds[ind]]]$filter_fn, filter_args)

            ## transform data back to original scale and store result in filtered_data
            if (!is.null(filter_control[[non_null_phi_inds[ind]]]$detransform_fn)) {
                filtered_data[, non_null_phi_inds[ind]] <-
                    filter_control[[non_null_phi_inds[ind]]]$detransform_fn(filtered_data_one_var)
            }
        }

        ## return
        return(cbind(data, as.data.frame(filtered_data)))
    }
}

#' Forward and reverse pass filter to correct phase shift introduced by one-pass
#' filter.  Uses stats::filter to do each filtering pass.  Closely based on
#' signal::filtfilt -- see the documentation there for further discussion.
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
two_pass_stats_filter <-
    function(filter, x, method, sides, circular) {
        #    y <- stats::filter(filter = filter, x = c(x, numeric(2 * length(filter))), method = method, sides = sides, circular = circular)
        y <-
            stats::filter(
                filter = filter,
                x = c(x, rep(x[length(x)], 2 * length(filter))),
                method = method,
                sides = sides,
                circular = circular
            )
        y <-
            rev(
                stats::filter(
                    filter = filter,
                    x = rev(y),
                    method = method,
                    sides = sides,
                    circular = circular
                )
            )[seq_along(x)]
        return(y)
    }

#' Forward and reverse pass filter to correct phase shift introduced by one-pass
#' filter.  Uses signal::filter to do each filtering pass.  Closely based on
#' signal::filtfilt -- see the documentation there for further discussion.  Two
#' main differences with that implementation:
#'   1) We impute internal missing values for the purposes of filtering, then
#'      replace them with NAs again after filtering
#'   2) We pad the end of the time series with the last observed value instead
#'      of zeros.  Based on very informal plotting, this seems to give better
#'      values at the end of the time series, which will be important for
#'      prediction
#'
#' @param filter vector of filter coefficients for a FIR filter
#' @param x vector of data to be filtered
#' @param method method for a call to stats::filter -- see the documentation at stats::filter
#' @param sides sides for a call to stats::filter -- see the documentation at stats::filter
#' @param circular circular for a call to stats::filter -- see the documentation at stats::filter
#'
two_pass_signal_filter <- function(filt, x, impute_fn) {
    ## Pad end with last observed value
    filt_Arma <- signal::as.Arma(filt)
    x_padded <- c(x,
                  rep(x[length(x)],
                      2 * max(
                          length(filt_Arma$a), length(filt_Arma$b)
                      )))

    ## Impute_fn should return a matrix with nrow = length(x_padded)
    ## and one column per imputation (for multiple imputation)
    x_imputed <- impute_fn(x_padded)

    ## Do two pass filtering on each colum of x_imputed
    y <- matrix(NA,
                nrow = nrow(x_imputed),
                ncol = ncol(x_imputed))
    for (j in seq_len(ncol(y))) {
        y[, j] <- signal::filter(filt = filt, x = x_imputed[, j])
        y[, j] <- rev(signal::filter(filt = filt, x = rev(y[, j])))
    }

    ## Take average of filtered values obtained with each imputed series
    ## Also, drop padded values at end
    y <- apply(y, 1, mean)[seq_along(x)]

    ## Re-insert NA values that were filled by imputation
    y[is.na(x)] <- NA

    ## Return
    return(y)
}

#' One-pass filter using signal::filter, with two additions:
#'   1) We impute internal missing values for the purposes of filtering, then
#'      replace them with NAs again after filtering
#'   2) We pad the end of the time series with the last observed value instead
#'      of zeros.  Based on very informal plotting, this seems to give better
#'      values at the end of the time series, which will be important for
#'      prediction
#'
#' @param filter vector of filter coefficients for a FIR filter
#' @param x vector of data to be filtered
#' @param method method for a call to stats::filter -- see the documentation at stats::filter
#' @param sides sides for a call to stats::filter -- see the documentation at stats::filter
#' @param circular circular for a call to stats::filter -- see the documentation at stats::filter
#'
one_pass_signal_filter <- function(filt, x, impute_fn) {
    ## Pad end with last observed value
    filt_Arma <- signal::as.Arma(filt)
    x_padded <- c(x,
                  rep(x[length(x)],
                      2 * max(
                          length(filt_Arma$a), length(filt_Arma$b)
                      )))

    ## Impute_fn should return a matrix with nrow = length(x_padded)
    ## and one column per imputation (for multiple imputation)
    x_imputed <- impute_fn(x_padded)

    ## Do two pass filtering on each colum of x_imputed
    y <- matrix(NA,
                nrow = nrow(x_imputed),
                ncol = ncol(x_imputed))
    for (j in seq_len(ncol(y))) {
        y[, j] <- signal::filter(filt = filt, x = x_imputed[, j])
    }

    ## Take average of filtered values obtained with each imputed series
    ## Also, drop padded values at end
    y <- apply(y, 1, mean)[seq_along(x)]

    ## Re-insert NA values that were filled by imputation
    y[is.na(x)] <- NA

    ## Return
    return(y)
}

#' Fill interior NAs in x by linear interpolation
#'
#' @param x a numeric vector of values to fill by linear interpolation
#'
#' @return a matrix with nrow = length(x) and one column, with interior NAs in x
#'   replaced by linearly interpolated values
interior_linear_interpolation <- function(x) {
    return(matrix(approx(seq_along(x), x, seq_along(x))$y))
}
