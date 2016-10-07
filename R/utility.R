## Miscellaneous utility functions
##
## update_vars_and_offsets (exported)
## compute_normalized_log_weights (exported)
## compute_offset_obs_vecs (exported)
## compute_na_rows_after_filter_and_offset (exported)
## log_score_loss (exported)
## neg_log_score_loss (exported)
## mase_from_kernel_weights_and_centers (exported)
## mae_from_kernel_weights_and_centers (exported)
## get_inds_smallest_k (exported)
## logspace_sub
## logspace_add
## logspace_sum
## logspace_sum_matrix_rows
## logspace_sub_matrix_rows
## mae (exported)
## mase (exported)

#' Update the list that keeps track of which combinations of variables and
#' offsets are included in the model by adding or removing (as appropriate) the
#' combination specified by update_var_name, update_offset_values, and
#' update_offset_type
#'
#' @param prev_vars_and_offsets list of previous variable/lag combinations to
#'     update
#' @param update_var_name name of the variable to update
#' @param update_offset_value value of the offset to update
#' @param update_offset_type type of the offset to update
#'
#' @return updated list of variable/lag combinations
#' @export
update_vars_and_offsets <- function(prev_vars_and_offsets,
                                    update_var_name,
                                    update_offset_value,
                                    update_offset_type) {
    updated_vars_and_offsets <- prev_vars_and_offsets

    existing_ind <-
        which(
            updated_vars_and_offsets$var_name == update_var_name &
                updated_vars_and_offsets$offset_value == update_offset_value &
                updated_vars_and_offsets$offset_type == update_offset_type
        )

    if (length(existing_ind) > 0) {
        ## remove variable/offset combination from model
        updated_vars_and_offsets <-
            updated_vars_and_offsets[-existing_ind, , drop = FALSE]
    } else {
        ## add variable/offset combination to model
        updated_vars_and_offsets <- rbind(
            updated_vars_and_offsets,
            data.frame(
                var_name = update_var_name,
                offset_value = as.integer(update_offset_value),
                offset_type = update_offset_type,
                combined_name = paste0(
                    update_var_name,
                    "_",
                    update_offset_type,
                    update_offset_value
                )
            )
        )

        ## Sort rows by combined_name to ensure uniqueness
        row_reordering <-
            order(updated_vars_and_offsets$combined_name)

        updated_vars_and_offsets <-
            updated_vars_and_offsets[row_reordering, , drop = FALSE]
    }

    ## Set rownames to NULL -- otherwise, the order in which variables are added
    ## and removed from the model shows up, and comparisons with identical() fail
    rownames(updated_vars_and_offsets) <- NULL

    return(updated_vars_and_offsets)
}

#' Normalize a vector of weights on log scale: given a vector of log_weights
#' such that exp(log_weights) are proportional to the final weights, update
#' so that exp(log_weights) sums to 1.
#'
#' @param log_weights: a vector of log(w) where w is proportional to the weights
#'
#' @return normalized log_weights so that sum(exp(log_weights)) = 1
#' @export
compute_normalized_log_weights <- function(log_weights) {
    ## normalize
    norm_const <-
        logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    log_weights <- log_weights - norm_const

    ## normalize again -- if the initial log_weights were "extreme", the norm_const
    ## computed above may be approximate.
    norm_const <-
        logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))

    return(log_weights - norm_const)
}

#' Compute a data frame with offset (lag and/or lead) observation vectors
#'
#' @param data a data frame
#' @param filter_control control for filtering data (compute_filter_values)
#' @param phi phi parameter for filtering (compute_filter_values)
#' @param vars_and_offsets a named list: The component name matches the name of
#'     one of the variables in data, and the component value is an integer
#'     vector of lags to include for that variable
#' @param time_name time names for adding to data columns
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be the largest lag used, but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#' @param trailing_rows_to_drop an integer specifying the number of trailing rows
#'     to drop.  This will typically be the largest prediction horizon used.
#' @param additional_rows_to_drop an integer vector specifying indices of
#'     additional rows to drop.  For example, if we are performing
#'     cross-validation, we might want to drop all rows within +/- 52 indices of
#'     the current prediction target.
#' @param na.action string specifying action to take with rows with NA. Supports
#'     "na.omit" (drops) and "na.pass" (do nothing)
#'
#' @return final data frame
#' @export
compute_offset_obs_vecs <- function(data,
                                    filter_control,
                                    phi,
                                    vars_and_offsets,
                                    time_name,
                                    leading_rows_to_drop = 0,
                                    trailing_rows_to_drop = 0,
                                    additional_rows_to_drop = NULL,
                                    na.action = "na.omit") {
    ## validate leading_rows_to_drop, trailing_rows_to_drop, additional_rows_to_drop
    if (any(
        c(
            leading_rows_to_drop,
            trailing_rows_to_drop,
            additional_rows_to_drop
        ) < 0 |
        c(
            leading_rows_to_drop,
            trailing_rows_to_drop,
            additional_rows_to_drop
        ) > nrow(data)
    )) {
        stop(
            "leading_rows_to_drop, trailing_rows_to_drop, and additional_rows_to_drop must be integers between 0 and nrow(data)"
        )
    }

    ## max_lag <- suppressWarnings(max(vars_and_offsets$offset_value[vars_and_offsets$offset_type == "lag"]))
    ## if (max_lag > -Inf && leading_rows_to_drop < max_lag) {
    ##     stop("leading_rows_to_drop must be at least as large as the largest lag in vars_and_offsets")
    ## }

    ## max_horizon <- suppressWarnings(max(vars_and_offsets$offset_value[vars_and_offsets$offset_type == "horizon"]))
    ## if (max_horizon > -Inf && trailing_rows_to_drop < max_horizon) {
    ##     stop("trailing_rows_to_drop must be at least as large as the largest horizon in vars_and_offsets")
    ## }

    ## rows to drop -- which rows should not be used as regression/density estimation examples
    ## either because the corresponding regression example cannot be formed
    ## or because we are performing cross-validation and don't want to use
    ## times adjacent to the prediction target
    if (is.logical(additional_rows_to_drop)) {
        additional_rows_to_drop <- which(additional_rows_to_drop)
    }
    rows_to_drop <-
        c(
            seq_len(leading_rows_to_drop),
            # leading -- too early for all lags to be available
            seq_len(trailing_rows_to_drop) + nrow(data) - trailing_rows_to_drop,
            # trailing -- too late for all prediction horizons to be available
            additional_rows_to_drop # additional (e.g., near prediction target)
        )

    ## create a data frame with one column for each entry in the vars_and_offsets argument
    result <- as.data.frame(matrix(
        NA,
        nrow = nrow(data),
        ncol = nrow(vars_and_offsets) +!(missing(time_name) ||
                                             is.null(time_name))
    ))
    if (missing(time_name) || is.null(time_name)) {
        colnames(result) <- vars_and_offsets$combined_name
    } else {
        colnames(result) <- c(vars_and_offsets$combined_name, time_name)
    }

    ## perform filtering if required
    if (!is.null(filter_control)) {
        filtered_data <- compute_filter_values(data = data,
                                               filter_control = filter_control,
                                               phi = phi)
    } else {
        filtered_data <- data
    }

    ## set column values in result
    for (new_var_ind in seq_len(nrow(vars_and_offsets))) {
        offset_val <- vars_and_offsets[new_var_ind, "offset_value"]
        offset_type <-
            as.character(vars_and_offsets[new_var_ind, "offset_type"])
        var_name <- vars_and_offsets[new_var_ind, "var_name"]
        combined_name <-
            as.character(vars_and_offsets[new_var_ind, "combined_name"])

        if (identical(offset_type, "lag")) {
            if (offset_val < nrow(result)) {
                result_inds <- seq(from = 1 + offset_val,
                                   to = nrow(result))
                data_inds <-
                    seq(from = 1,
                        to = nrow(result) - offset_val)
            } else {
                result_inds <- c()
                data_inds <- c()
            }
        } else {
            # assumed only alternative is "horizon"
            if (offset_val < nrow(result)) {
                result_inds <- seq(from = 1,
                                   to = nrow(result) - offset_val)
                data_inds <-
                    seq(from = 1 + offset_val,
                        to = nrow(result))
            } else {
                result_inds <- c()
                data_inds <- c()
            }
        }

        result[result_inds, combined_name] <-
            filtered_data[data_inds, var_name]
    }

    if (!(missing(time_name) || is.null(time_name))) {
        result[, time_name] <- data[, time_name]
    }

    ## drop specified rows
    if (length(rows_to_drop) > 0) {
        result <- result[-rows_to_drop, , drop = FALSE]
    }

    ## Check NAs
    if (identical(na.action, "na.omit")) {
        na_rows <- which(apply(result, 1, anyNA))
        if (length(na_rows) > 0) {
            result <- result[-na_rows, , drop = FALSE]
            if (nrow(result) == 0) {
                stop(
                    "0 rows in cross validation data after filtering, computing offsets and dropping NA rows"
                )
            }
        }
    } else if (!identical(na.action, "na.pass")) {
        stop("Unsupported na.action")
    }

    return(result)
}

#' Compute rows of the data frame that contain NA values after performing the
#' filtering and observation offset steps.  Also include values of
#' prediction_inds_not_included that were specified in kcde_control, if
#' applicable.
#'
#' @param data data frame of observations to filter and offset
#' @param phi list of parameters for filtering
#' @param vars_and_offsets data frame describing variables to offset and size
#'     and direction of offset
#' @param kcde_control list of control parameters for kcde
#'
#' @return integer vector with rows of data that contain NA values after
#'     filtering and offsetting observations
#' @export
compute_na_rows_after_filter_and_offset <- function(data,
                                                    phi,
                                                    vars_and_offsets,
                                                    kcde_control) {
    phi_init <- initialize_phi(
        prev_phi = phi,
        updated_vars_and_offsets = vars_and_offsets,
        update_var_name = vars_and_offsets$var_name,
        update_offset_value = vars_and_offsets$offset_value,
        update_offset_type = vars_and_offsets$offset_type,
        data = data,
        kcde_control = kcde_control
    )

    filtered_and_lagged_data <- compute_offset_obs_vecs(
        data = data,
        filter_control = kcde_control$filter_control,
        phi = phi_init,
        vars_and_offsets = vars_and_offsets,
        time_name = kcde_control$time_name,
        leading_rows_to_drop = 0L,
        trailing_rows_to_drop = 0L,
        additional_rows_to_drop = NULL,
        na.action = "na.pass"
    )

    all_na_drop_rows <- which(apply(filtered_and_lagged_data[vars_and_offsets$combined_name],
                                    1,
                                    anyNA))

    if (!is.null(kcde_control$prediction_inds_not_included)) {
        all_na_drop_rows <- unique(c(
            all_na_drop_rows,
            kcde_control$prediction_inds_not_included
        ))
    }

    return(all_na_drop_rows)
}


### Loss Functions

#' Compute log score based on a vector of log(p(data)) where p is a predictive
#' distribution
#' @export
log_score_loss <- function(prediction_result, ...) {
    return(sum(prediction_result))
}

#' Compute negative log score based on a vector of log(p(data)) where p is a
#' predictive distribution
#' @export
neg_log_score_loss <- function(prediction_result, ...) {
    return(-1 * sum(prediction_result))
}

#' Compute MASE from predictions that are in the form of kernel weights and
#' centers.  This is currently broken because we need predictions for a whole
#' time series to compute MASE, and the current "implementation" only takes
#' values for one time point.
## mase_from_kernel_weights_and_centers <- function(kernel_weights_and_centers,
##                                                  obs) {
##     stop("mase_from_kernel_weights_and_centers is broken -- do not use")
##     pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
##     mase(obs, pred)
## }

#' Compute "mean" absolute error for one time point from prediction in the form
#' of kernel weights and centers.
#'
#' @param kernel_weights_and_centers a named list with two components: weights
#'     is a vector of kernel weights, and centers is a vector of kernel centers
#' @param obs is a numeric with length 1 containing the observed value for one
#'     time point
#'
#' @return abs(obs - prediction) where prediction is the weighted mean of the
#'     kernel centers.
#' @export
mae_from_kernel_weights_and_centers <-
    function(kernel_weights_and_centers,
             obs) {
        pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
        mae(obs, pred)
    }

#' Get the indices of the smallest k elements of v.
#'
#' @param v a vector
#' @param k number of indices to return
#'
#' @return a vector of length k containing the indices of the k smallest
#'     elements of v, in ascending order.
#' @export
get_inds_smallest_k <- function(v, k) {
    take <- if (k > length(v))
        length(v)
    else
        k
    return(order(v, decreasing = FALSE)[seq_len(take)])
}


### Numerical Functions

## Interface to R's C API for logspace arithmetic

#' Subtract two values in logspace
#'
#' @param logx minuend
#' @param logy subtrahend
#'
#' @return result of the subtraction
#' @useDynLib kcde logspace_sub_C
logspace_sub <- function(logx, logy) {
    return(.Call("logspace_sub_C", as.numeric(logx), as.numeric(logy)))
}

#' Add two values in logspace
#'
#' @param logx first number
#' @param logy second number
#'
#' @return sum of the numbers
#' @useDynLib kcde logspace_add_C
logspace_add <- function(logx, logy) {
    return(.Call("logspace_add_C", as.numeric(logx), as.numeric(logy)))
}

#' Find sum of numeric vector in logspace
#'
#' @param logx the vector to sum
#'
#' @return sum of all the elements in the vector
logspace_sum <- function(logx) {
    dim(logx) <- c(1, length(logx))
    return(logspace_sum_matrix_rows(logx))
}

#' Row-wise sum of the given matrix in logspace
#'
#' @param logX the matrix to sum
#'
#' @return a vector with each element summing corresponding row of matrix
#' @useDynLib kcde logspace_sum_matrix_rows_C
logspace_sum_matrix_rows <- function(logX) {
    return(.Call(
        "logspace_sum_matrix_rows_C",
        as.numeric(logX),
        as.integer(nrow(logX)),
        as.integer(ncol(logX))
    ))
}

#' Row-wise difference of matrix with two columns in logspace
#'
#' @param logX the matrix
#'
#' @return a vector with each element being the difference of corresponding
#'     elements in the first and second column
#' @useDynLib kcde logspace_sub_matrix_rows_C
logspace_sub_matrix_rows <- function(logX) {
    if (!is.matrix(logX) || !identical(ncol(logX), 2L))
        stop("logX must be a matrix with 2 columns")

    return(.Call(
        "logspace_sub_matrix_rows_C",
        as.numeric(logX),
        as.integer(nrow(logX))
    ))
}


### Forecast Evaluation

#' Compute Mean Absolute Error to evaluate point predictions
#'
#' @param obs vector of observed values
#' @param pred vector of point predictions
#'
#' @return mean absolute error for given vectors
#' @export
mae <- function(obs, pred) {
    if (length(obs) != length(pred)) {
        stop("arguments must have the same length")
    }
    return(mean(abs(obs - pred)))
}

#' Compute Mean Absolute Scaled Error (MASE) to evaluate point predictions
#'
#' @param obs vector of observed values
#' @param pred vector of point predictions
#'
#' @return mean absolute scaled error for given vectors
#' @export
mase <- function(obs, pred) {
    mean_forecast_error <- mae(obs, pred)
    mean_naive_error <- mae(obs[-length(obs)], obs[-1])
    return(mean_forecast_error / mean_naive_error)
}
