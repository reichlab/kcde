## Functions to set up the estimation routine for kcde
##
## create_kcde_control
## create_kcde_control_default
## get_default_kernel_components
## validate_kcde_control


#' Assemble a list of kcde_control parameters for the kcde function with
#'     user-specified values.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param prediction_horizons integer vector: the number of time steps between
#'     the last observation and the time at which we make a prediction
#' @param kernel_components a list with one component for each component of the
#'     kernel function.  Each component is a list with the following entries:
#'       - vars_and_lags: a data frame with two columns: "var_name" and
#'         "lag_value".  Each row specifies a combination of variable and lag
#'         that is included in this component of the kernel function.
#'       - kernel_fn: a function to evaluate the kernel
#'       - theta_fixed: a named list of parameters to kernel_fn whose values are
#'         held fixed (i.e., not estimated)
#'       - theta_est: a named list of parameters to kernel_fn whose values are
#'         to be estimated
#'       - initialize_theta_fn: a function to initialize both theta_fixed and
#'         theta_est
#'       - initialize_theta_args: a named list of arguments to
#'         initialize_theta_fn
#'       - vectorize_theta_est: a function that converts theta_est into an
#'         ordered vector on a scale suitable for passing as the first argument
#'         to optim.  Required to return a list with three components:
#'           (1) theta_est - vector of parameters to be estimated
#'           (2) lb - vector of lower bounds to theta_est
#'           (3) ub - vector of upper bounds to theta_est
#'       - update_theta_from_vectorized_theta_est: a function that updates
#'         theta_est (in list form) from theta_est (in vector form).
#' @param crossval_buffer during cross-validation, the number of indices before
#'     the time at which we are making a prediction to drop from the "training
#'     examples".
#' @param loss_fn_name a string giving the name of the function use to
#'     compute loss from predictions
#' @param loss_fn_args a named list giving arguments to the loss function
#' @param par_packages a character vector containing names of packages that need
#'     to be loaded in instances of R when computations are performed in
#'     parallel.
#' 
#' @return the (at this point, unvalidated) list of kcde_control parameters
create_kcde_control <- function(X_names,
        y_names,
        time_name,
        prediction_horizons,
        kernel_components,
        filter_control,
        crossval_buffer,
        loss_fn,
        loss_fn_prediction_args,
        loss_args,
        par_packages = NULL,
        na.action = "na.omit") {
    kcde_control <- list()
    
    kcde_control$X_names <- X_names
    kcde_control$y_names <- y_names
    kcde_control$time_name <- time_name
    
    kcde_control$prediction_horizons <- prediction_horizons
    
    kcde_control$kernel_components <- kernel_components
    kcde_control$filter_control <- filter_control
    
    kcde_control$crossval_buffer <- crossval_buffer
    
    kcde_control$loss_fn <- loss_fn
    kcde_control$loss_fn_prediction_args <- loss_fn_prediction_args
    kcde_control$loss_args <- loss_args
    
    kcde_control$par_packages <- par_packages
    
    kcde_control$na.action <- na.action
    
    return(kcde_control)
}

#' Assemble a list of kcde_control parameters for the kcde function with default
#'     values
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param data a data frame where rows are consecutive observations
#' 
#' @return the list of kcde_control parameters
create_kcde_control_default <- function(X_names, y_names, data) {
    kcde_control <- list()
    
    kcde_control$X_names <- X_names
    kcde_control$y_names <- y_names
    
    kcde_control$kernel_components <- get_default_kernel_components(X_names,
        y_names,
        time_name,
        data)
    
    kcde_control$loss_fn_name <- "mase"
    kcde_control$loss_fn_args <- list()
    
    return(kcde_control)
}

#' Get default kernel functions based on a brief look at the data.  This is
#' unreliable.  Update to return periodic_kernel if X_names[i] == time_name?
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return a list of default parameters for kernel components -- probably all bad
get_default_kernel_components <- function(X_names, y_names, data) {
	stop("Function get_default_kernel_components is not yet implemented")

    return(kernel_components)
}

#' Validate kcde_control parameters for kcde -- not implemented
#' 
#' @param kcde_control a list of kcde_control parameters for kcde
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return no return value -- either stops with an error or not.
validate_kcde_control <- function(kcde_control, X_names, y_names, time_name, data) {
#    warning("kcde kcde_control parameter validation not yet implemented")
}
