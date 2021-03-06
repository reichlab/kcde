## Functions to perform prediction from new data and an object with a fit
##
## kcde_predict
## kcde_predict_given_lagged_obs
## kcde_kernel_centers_and_weights_predict_given_lagged_obs
## kcde_point_predict_given_lagged_obs
## kcde_dist_predict_given_lagged_lead_obs


#' Make predictions from an estimated kcde model forward prediction_horizon time
#' steps from the end of predict_data, based on the weighting variables, lags,
#' kernel functions, and bandwidths specified in the kcde_fit object.
#' 
#' @param kcde_fit is an object representing a fitted kcde model
#' @param prediction_data is a vector of data points to use in prediction
#' @param normalize_weights boolean, should the weights be normalized?
#' @param prediction_type character; either "distribution" or "point",
#'     indicating the type of prediction to perform.
#' 
#' @return an object with prediction results; the contents depend on the value
#'     of prediction_type
kcde_predict <- function(kcde_fit,
        prediction_data,
        leading_rows_to_drop = max(kcde_fit$vars_and_offsets$offset_value[kcde_fit$vars_and_offsets$offset_type == "lag"]),
        trailing_rows_to_drop = max(kcde_fit$vars_and_offsets$offset_value[kcde_fit$vars_and_offsets$offset_type == "horizon"]),
        additional_training_rows_to_drop = NULL,
        prediction_type = "distribution",
        n,
        p,
        q,
        prediction_test_lead_obs,
        log = FALSE) {
    ## get training and prediction examples
    ## create data frame of "examples" -- lagged observation vectors and
    ## corresponding prediction targets
    training_examples <- compute_offset_obs_vecs(data = kcde_fit$train_data,
        filter_control = kcde_fit$kcde_control$filter_control,
        phi = kcde_fit$phi_hat,
        vars_and_offsets = kcde_fit$vars_and_offsets,
        time_name = kcde_fit$kcde_control$time_name,
        leading_rows_to_drop = leading_rows_to_drop,
        trailing_rows_to_drop = trailing_rows_to_drop,
        additional_rows_to_drop = kcde_fit$kcde_control$prediction_inds_not_included,
        na.action = kcde_fit$kcde_control$na.action)
    
    ## Keep only non-NA rows in training data
    training_non_na_rows <- which(!apply(
            training_examples[kcde_fit$vars_and_offsets$combined_name],
            1,
            anyNA))
    
    training_examples <- training_examples[training_non_na_rows, , drop = FALSE]
    
    prediction_examples <- compute_offset_obs_vecs(data = prediction_data,
        filter_control = kcde_fit$kcde_control$filter_control,
        phi = kcde_fit$phi_hat,
        vars_and_offsets = kcde_fit$vars_and_offsets[kcde_fit$vars_and_offsets$offset_type == "lag", , drop = FALSE],
        time_name = kcde_fit$kcde_control$time_name,
        leading_rows_to_drop = 0L,
        trailing_rows_to_drop = 0L,
        additional_rows_to_drop = NULL,
        na.action = "na.pass")
    
    ## Keep only last row in prediction data
    prediction_examples <- prediction_examples[nrow(prediction_examples), , drop = FALSE]
    
    predictive_var_combined_names <- kcde_fit$vars_and_offsets$combined_name[kcde_fit$vars_and_offsets$offset_type == "lag"]
    target_var_combined_names <- kcde_fit$vars_and_offsets$combined_name[kcde_fit$vars_and_offsets$offset_type == "horizon"]
    
    ## At present, can't handle missing values in predictive variables
    if(any(is.na(prediction_examples[, predictive_var_combined_names]))) {
      warning("Can't predict based on missing values!")
      return(NA)
    }
    
    ## do prediction
    return(kcde_predict_given_lagged_obs(
            train_lagged_obs = training_examples[, predictive_var_combined_names, drop = FALSE],
            train_lead_obs = training_examples[, target_var_combined_names, drop = FALSE],
            prediction_lagged_obs = prediction_examples[, predictive_var_combined_names, drop = FALSE],
#            prediction_test_lead_obs = prediction_examples[, target_var_combined_names, drop = FALSE],
            prediction_test_lead_obs = prediction_test_lead_obs,
            kcde_fit = kcde_fit,
            prediction_type = prediction_type,
            n = n,
            p = p,
            q = q,
            log = log))
}

#' Make predictions from an estimated kcde model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the kcde_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' @param prediction_type character; either "distribution" or "point",
#'     indicating the type of prediction to perform.
#' 
#' @return an object with prediction results; the contents depend on the value
#'     of prediction_type
kcde_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    prediction_test_lead_obs,
    kcde_fit,
    prediction_type="distribution",
    n,
    p,
    q,
    log = FALSE) {
    
    if(identical(prediction_type, "centers-and-weights")) {
        kernel_centers_and_weights <-
            kcde_kernel_centers_and_weights_predict_given_lagged_obs(
                train_lagged_obs = train_lagged_obs,
                train_lead_obs = train_lead_obs,
                prediction_lagged_obs = prediction_lagged_obs,
                kcde_fit = kcde_fit)
        
        return(kernel_centers_and_weights)
    } else if(identical(prediction_type, "distribution")) {
    	return(kcde_dist_predict_given_lagged_lead_obs(train_lagged_obs,
            train_lead_obs,
            prediction_lagged_obs,
            prediction_test_lead_obs,
            kcde_fit,
            log = log))
    } else if(identical(prediction_type, "point")) {
        kernel_centers_and_weights <-
            kcde_kernel_centers_and_weights_predict_given_lagged_obs(
                train_lagged_obs = train_lagged_obs,
                train_lead_obs = train_lead_obs,
                prediction_lagged_obs = prediction_lagged_obs,
                kcde_fit = kcde_fit)
        
        return(kcde_point_predict_given_kernel_centers_and_weights(
            kernel_centers_and_weights = kernel_centers_and_weights,
            kcde_fit = kcde_fit))
    } else if(identical(prediction_type, "sample")) {
        return(kcde_sample_predict_given_lagged_lead_obs(
            n = n,
            train_lagged_obs = train_lagged_obs,
            train_lead_obs = train_lead_obs,
            prediction_lagged_obs = prediction_lagged_obs,
            kcde_fit = kcde_fit))
    } else if(identical(prediction_type, "quantile")) {
        return(kcde_quantile_predict_given_lagged_lead_obs(
                p = p,
                n = n,
                train_lagged_obs = train_lagged_obs,
                train_lead_obs = train_lead_obs,
                prediction_lagged_obs = prediction_lagged_obs,
                kcde_fit = kcde_fit))
    } else if(identical(prediction_type, "prob")) {
        return(kcde_prob_predict_given_lagged_lead_obs(
                q = q,
                n = n,
                train_lagged_obs = train_lagged_obs,
                train_lead_obs = train_lead_obs,
                prediction_lagged_obs = prediction_lagged_obs,
                kcde_fit = kcde_fit))
    } else {
        stop("Invalid prediction type.")
    }
}

#' Make predictions from an estimated kcde model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the kcde_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a named list with four components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
#'     pred_bws: a named vector of bandwidths, one per column of centers
kcde_kernel_centers_and_weights_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    kcde_fit) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    unnormalized_log_weights <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        kernel_components = kcde_fit$kcde_control$kernel_components,
        theta = kcde_fit$theta_hat,
        log = TRUE)
    
    log_weights <- compute_normalized_log_weights(unnormalized_log_weights)
    
    return(list(unnormalized_log_weights=unnormalized_log_weights,
        log_weights=log_weights,
        weights=exp(log_weights),
        conditioning_vars=train_lagged_obs,
        centers=train_lead_obs))
}


#' Make predictions from an estimated kcde model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the kcde_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     one row, representing one time point.  Each column is a (lagged) variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
kcde_dist_predict_given_lagged_lead_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    prediction_test_lead_obs,
    kcde_fit,
    log) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    log_kernel_values_x <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        kernel_components = kcde_fit$kcde_control$kernel_components,
        theta = kcde_fit$theta_hat,
        log = TRUE)
    
    ## compute log(sum(kernel values x))
    log_sum_kernel_values_x <- logspace_sum(log_kernel_values_x)
    
    log_result <- sapply(seq_len(nrow(prediction_test_lead_obs)),
        function(prediction_test_row_ind) {
            log_kernel_values_xy <- compute_kernel_values(
                cbind(train_lagged_obs, train_lead_obs),
                cbind(prediction_lagged_obs,
                    prediction_test_lead_obs[prediction_test_row_ind, , drop = FALSE]),
                kernel_components = kcde_fit$kcde_control$kernel_components,
                theta = kcde_fit$theta_hat,
                log = TRUE
            )
            
            return(logspace_sum(log_kernel_values_xy) -
                    log_sum_kernel_values_x)
        })
    
    if(log) {
        return(log_result)
    } else {
        return(exp(log_result))
    }
}

#' Make predictions from an estimated kcde model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the kcde_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param prediction_test_lead_obs is a matrix (with column names) containing
#'     prediction target vectors computed from the prediction data.  Each row
#'     represents one time point.  Each column is a (leading) target variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
kcde_point_predict_given_kernel_centers_and_weights <- function(kernel_centers_and_weights,
    kcde_fit) {
    stop("point predictions for kcde are not yet implemented")
}

#' Draw a sample from the predictive distribution corresponding to an estimated
#' kcde model forward prediction_horizon time steps from the end of predict_data.
#' This function requires that the kernel weights and centers have already been
#' computed.
#' 
#' @param n sample size for sample from predictive distribution used in
#'     approximating quantiles
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' 
#' @return a matrix with samples from the predictive distribution
kcde_sample_predict_given_lagged_lead_obs <- function(n,
    train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    kcde_fit) {
    
    kernel_centers_and_weights <-
        kcde_kernel_centers_and_weights_predict_given_lagged_obs(
            train_lagged_obs = train_lagged_obs,
            train_lead_obs = train_lead_obs,
            prediction_lagged_obs = prediction_lagged_obs,
            kcde_fit = kcde_fit)
    
    result <- matrix(NA, nrow = n, ncol = 1)
    sampled_kernel_inds <- sample(length(kernel_centers_and_weights$weights),
        size = n,
        replace = TRUE,
        prob = kernel_centers_and_weights$weights)
    
    for(kernel_ind in unique(sampled_kernel_inds)) {
        result_inds <- which(sampled_kernel_inds == kernel_ind)
        
        result[result_inds, ] <- simulate_values_from_product_kernel(n = length(result_inds),
            conditioning_obs = prediction_lagged_obs,
            center = cbind(kernel_centers_and_weights$conditioning_vars[kernel_ind, , drop = FALSE],
                kernel_centers_and_weights$centers[kernel_ind, , drop = FALSE]),
            kernel_components = kcde_fit$kcde_control$kernel_components,
            theta = kcde_fit$theta)[, colnames(kernel_centers_and_weights$centers)]
    }
    
    return(result)
}

#' Draw a sample from the predictive distribution corresponding to an estimated
#' kcde model forward prediction_horizon time steps from the end of predict_data.
#' This function requires that the kernel weights and centers have already been
#' computed.
#' 
#' @param p vector of probabilities which to compute quantiles
#' @param n sample size for sample from predictive distribution used in
#'     approximating quantiles
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' 
#' @return a vector of quantiles
kcde_quantile_predict_given_lagged_lead_obs <- function(
    p,
    n,
    train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    kcde_fit) {
    
    predictive_sample <- kcde_sample_predict_given_lagged_lead_obs(
        n = n,
        train_lagged_obs = train_lagged_obs,
        train_lead_obs = train_lead_obs,
        prediction_lagged_obs = prediction_lagged_obs,
        kcde_fit = kcde_fit)
    
    return(apply(predictive_sample, 2, quantile, probs = p))
}

#' Compute the conditional (predictive) probability that the predictive target variables
#' are all less than or equal to the corresponding elements of the rows of q.
#' 
#' @param p vector of probabilities which to compute quantiles
#' @param n sample size for sample from predictive distribution used in
#'     approximating quantiles
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param kcde_fit is an object representing a fitted kcde model
#' 
#' @return a vector of probabilities
kcde_prob_predict_given_lagged_lead_obs <- function(
    q,
    n,
    train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    kcde_fit) {
    
    predictive_sample <- kcde_sample_predict_given_lagged_lead_obs(
        n = n,
        train_lagged_obs = train_lagged_obs,
        train_lead_obs = train_lead_obs,
        prediction_lagged_obs = prediction_lagged_obs,
        kcde_fit = kcde_fit)
    
    if(!is.matrix(q)) {
        q <- as.matrix(q)
    }
    
    return(apply(q, 1, function(q_row) {
        inds <- apply(predictive_sample, 1, function(sample_row) {
            return(all(sample_row < q_row))
        })
        return(sum(inds) / nrow(predictive_sample))
    }))
}
