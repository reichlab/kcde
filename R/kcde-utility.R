## Miscellaneous utility functions
##
## update_vars_and_lags
## assemble_training_examples
## assemble_prediction_examples
## compute_normalized_log_weights
## compute_lagged_obs_vecs
## mase_from_kernel_weights_and_centers
## mae_from_kernel_weights_and_centers
## get_inds_smallest_k
## logspace_sub
## logspace_add
## logspace_sum
## logspace_sum_matrix_rows
## logspace_sub_matrix_rows
## mase
## mae

#' Update the list that keeps track of which combinations of variables and lags
#' are included in the model by adding or removing (as appropriate) the
#' combination specified by update_var_name and update_lag_values
#'
#' @param prev_vars_and_lags list of previous variable/lag combinations to update
#' @param update_var_name name of the variable to update
#' @param update_lag_value value of the lag to update
#'
#' @return updated list of variable/lag combinations
update_vars_and_lags <- function(prev_vars_and_lags,
		update_var_name,
		update_lag_value) {
	updated_vars_and_lags <- prev_vars_and_lags
	
	existing_ind <- which(prev_vars_and_lags$var_name == update_var_name &
		prev_vars_and_lags$lag_value == update_lag_value)
	
    if(length(existing_ind) > 0) {
        ## remove variable/lag combination from model
		updated_vars_and_lags <- updated_vars_and_lags[-existing_ind, , drop = FALSE]
    } else {
        ## add variable/lag combination to model
		updated_vars_and_lags <- rbind(updated_vars_and_lags,
			data.frame(var_name = update_var_name,
				lag_value = as.integer(update_lag_value),
				combined_name = paste0(update_var_name,
					"_lag",
					update_lag_value)))
    }
    
    return(updated_vars_and_lags)
}

#' Construct data frame of lagged observation vectors and data frame of
#' corresponding prediction targets.
#' 
#' @param data data frame with observations of all variables at all times
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param y_names names of variables included as prediction targets
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#' @param additional_rows_to_drop an integer vector specifying indices of
#'     additional rows to drop.  For example, if we are performing
#'     cross-validation, we might want to drop all rows within +/- 52 indices of
#'     the current prediction target.
#' @param prediction_horizons an integer vector specifying the number of time
#'     steps between the last observation and the prediction target.
#' @param drop_trailing_rows boolean:  drop the last prediction_horizon rows?
#'     These are the rows for which we can form lagged observation vectors, but
#'     we cannot obtain a corresponding prediction target.
#'     
#' @return a list with two components: lagged_obs is a data frame with lagged
#'     observation vectors, and lead_obs is a data frame with corresponding
#'     prediction targets.
assemble_training_examples <- function(data,
		vars_and_lags,
    	y_names,
    	leading_rows_to_drop,
    	additional_rows_to_drop,
    	prediction_horizons,
    	drop_trailing_rows=TRUE) {
    ## which rows should not be used as regression/density estimation examples
    ## either because the corresponding regression example cannot be formed
    ## or because we are performing cross-validation and don't want to use
    ## times adjacent to the prediction target

    ## too early
    all_train_rows_to_drop <- seq_len(leading_rows_to_drop)

    ## passed in indices -- near prediction target
    all_train_rows_to_drop <- c(all_train_rows_to_drop, additional_rows_to_drop)
    
    ## too late
    max_prediction_horizon <- max(prediction_horizons)
    if(drop_trailing_rows) {
        all_train_rows_to_drop <- c(all_train_rows_to_drop,
            seq(from=nrow(data) - max_prediction_horizon + 1,
                to=nrow(data)))
    }
    
    ## compute lagged observation vectors for train and prediction data
    train_lagged_obs <- compute_lagged_obs_vecs(data,
		vars_and_lags,
        all_train_rows_to_drop)
    
    ## compute lead observation series for train data
    train_lead_obs <- data.frame(rep(NA, nrow(data)))
    for(y_name in y_names) {
        for(prediction_horizon in prediction_horizons) {
            train_lead_obs_inds <-
                seq(from=1, to=nrow(data)) + prediction_horizon
            component_name <- paste0(y_name, "_horizon", prediction_horizon)
            train_lead_obs[[component_name]] <-
                data[train_lead_obs_inds, y_name, drop=TRUE]
        }
    }
    train_lead_obs <- train_lead_obs[-all_train_rows_to_drop, , drop=FALSE]
    
    ## drop initial column of NA's
    train_lead_obs <- train_lead_obs[, -1, drop=FALSE]
    
    return(list(lagged_obs=train_lagged_obs,
        lead_obs=train_lead_obs))
}

#' Construct data frame of lagged observation vectors.
#' 
#' @param data data frame with observations of all variables at all times
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#'     
#' @return a list with one component: lagged_obs is a data frame with lagged
#'     observation vectors.
assemble_prediction_examples <- function(data,
		vars_and_lags,
    	y_names,
    	prediction_horizons,
    	leading_rows_to_drop) {
    all_prediction_rows_to_drop <- seq_len(leading_rows_to_drop) # too early

    prediction_lagged_obs <- compute_lagged_obs_vecs(data,
		vars_and_lags,
        all_prediction_rows_to_drop)
    
    prediction_lead_obs <- data.frame(rep(NA, nrow(data)))
    for(y_name in y_names) {
        for(prediction_horizon in prediction_horizons) {
            prediction_lead_obs_inds <-
                seq(from=1, to=nrow(data)) + prediction_horizon
            component_name <- paste0(y_name, "_horizon", prediction_horizon)
            prediction_lead_obs[[component_name]] <-
                data[prediction_lead_obs_inds, y_name, drop=TRUE]
        }
    }
    prediction_lead_obs <- prediction_lead_obs[-all_prediction_rows_to_drop, , drop=FALSE]
    
    return(list(lagged_obs=prediction_lagged_obs))
}

#' Normalize a vector of weights on log scale: given a vector of log_weights
#' such that exp(log_weights) are proportional to the final weights, update
#' so that exp(log_weights) sums to 1.
#' 
#' @param log_weights: a vector of log(w) where w is proportional to the weights
#' 
#' @return normalized log_weights so that sum(exp(log_weights)) = 1
compute_normalized_log_weights <- function(log_weights) {
    ## normalize
    norm_const <- logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    log_weights <- log_weights - norm_const
    
    ## normalize again -- if the initial log_weights were "extreme", the norm_const
    ## computed above may be approximate.
    norm_const <- logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    
    return(log_weights - norm_const)
}


#' Compute a data frame with lagged observation vectors
#' 
#' @param data a data frame
#' @param vars_and_lags a named list: The component name matches the name of 
#'     one of the variables in data, and the component value is an integer
#'     vector of lags to include for that variable
#' @param rows_to_drop an integer vector specifying rows to drop after computing
#'     lagged observation vectors.
compute_lagged_obs_vecs <- function(data,
		vars_and_lags,
    	rows_to_drop) {
    ## set or validate leading_rows_to_drop
    max_lag <- max(vars_and_lags$lag_value)
    if(any(rows_to_drop < 0 | rows_to_drop > nrow(data))) {
        stop("all entries of rows_to_drop must integers between 1 and nrow(data)")
    } else if(!all(seq_len(max_lag) %in% rows_to_drop)) {
        stop("all integers between 1 and the maximum entry in the lags argument must be contained in rows_to_drop")
    }
    
    ## create a data frame with one column for each entry in the lags argument
    result <- as.data.frame(matrix(NA,
        nrow=nrow(data),
        ncol=nrow(vars_and_lags)))
    colnames(result) <- vars_and_lags$combined_name
    
    ## set column values in result
	for(new_var_ind in seq_len(nrow(vars_and_lags))) {
		lag_val <- vars_and_lags[new_var_ind, "lag_value"]
		var_name <- as.character(vars_and_lags[new_var_ind, "var_name"])
		combined_name <- as.character(vars_and_lags[new_var_ind, "combined_name"])
		
		result_inds <- seq(from=lag_val + 1, to=nrow(result))
		data_inds <- seq(from=1, to=nrow(result) - lag_val)
		
		result[result_inds, combined_name] <- data[data_inds, var_name]
	}
	
    ## drop specified rows
    if(length(rows_to_drop) > 0) {
        result <- result[-rows_to_drop, , drop=FALSE]
    }
    
    return(result)
}




### loss functions

#' Compute MASE from predictions that are in the form of kernel weights and
#' centers.  This is currently broken because we need predictions for a whole
#' time series to compute MASE, and the current "implementation" only takes
#' values for one time point....
mase_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    stop("mase_from_kernel_weights_and_centers is broken -- do not use")
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mase(obs, pred)
}

#' Compute "mean" absolute error for one time point from prediction in the form
#' of kernel weights and centers.
#' 
#' @param kernel_weights_and_centers a named list with two components: weights
#'     is a vector of kernel weights, and centers is a vector of kernel centers
#' @param obs is a numeric with length 1 containign the observed value for one
#'     time point
#' 
#' @return abs(obs - prediction) where prediction is the weighted mean of the
#'     kernel centers.
mae_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mae(obs, pred)
}

#' Get the indices of the smallest k elements of v.  This code currently assumes
#' that k >= length(v)
#' 
#' @param v a vector
#' @param k number of indices to return
#' 
#' @return a vector of length k containing the indices of the k smallest
#'   elements of v, in ascending order.
get_inds_smallest_k <- function(v, k) {
    return(order(v, decreasing=FALSE)[seq_len(k)])
}



### numerical functions

## interface to R's C API for logspace arithmetic

logspace_sub <- function(logx, logy) {
    return(.Call("logspace_sub_C", as.numeric(logx), as.numeric(logy)))
}

logspace_add <- function(logx, logy) {
    return(.Call("logspace_add_C", as.numeric(logx), as.numeric(logy)))
}

logspace_sum <- function(logx) {
    dim(logx) <- c(1, length(logx))
    return(logspace_sum_matrix_rows(logx))
}

logspace_sum_matrix_rows <- function(logX) {
    return(.Call("logspace_sum_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX)), as.integer(ncol(logX))))
}

logspace_sub_matrix_rows <- function(logX) {
    if(!is.matrix(logX) || !identical(ncol(logX), 2L))
        stop("logX must be a matrix with 2 columns")
    
    return(.Call("logspace_sub_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX))))
}


### forecast evaluation

#' Compute Mean Absolute Scaled Error (MASE) to evaluate point predictions
#' 
#' @param obs vector of observed values
#' @param pred vector of point predictions
mase <- function(obs, pred) {
    mean_forecast_error <- mean(abs(obs - pred))
    mean_naive_error <- mean(abs(obs[-length(obs)] - obs[-1]))
    return(mean_forecast_error / mean_naive_error)
}

#' Compute Mean Absolute Error to evaluate point predictions
#' 
#' @param obs vector of observed values
#' @param pred vector of point predictions
mae <- function(obs, pred) {
    if(length(obs) != length(pred)) {
        stop("obs and pred must have the same length")
    }
    
    return(mean(abs(obs - pred)))
}
