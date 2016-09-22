## pdtmvn kernel parameterized by mode and variance function
## 
## log_pdtmvn_mode_centered_kernel
## rlog_pdtmvn_mode_centered_kernel
## vectorize_params_log_pdtmvn_mode_centered_kernel
## update_theta_from_vectorized_theta_est_log_pdtmvn_mode_centered_kernel
## initialize_params_log_pdtmvn_mode_centered_kernel

#' Evaluate the kernel function given by the log-pdtmvn mode centered distribution.
#' 
#' @param a matrix of values at which to evaluate the kernel function, with
#'     column names specified.  Each row is an observation, each column is an
#'     observed variable.
#' @param center a real vector, center point for the kernel function
#' @param bw bandwidth matrix
#' @param bw_continuous the portion of bw corresponding to continuous variables
#' @param conditional_bw_discrete the Schur complement of the portion of bw
#'     corresponding to discrete variables
#' @param conditional_center_discrete_offset_multiplier Sigma_dc Sigma_c^{-1}.
#'     This is used in computing the mean of the underlying multivariate normal
#'     distribution for the discrete variables conditioning on the continuous
#'     variables.
#' @param continuous_vars character vector with names of variables that are to
#'     be treated as continuous.  May contain entries that do not appear in
#'     colnames(x).
#' @param discrete_vars character vector with names of variables that are to be
#'     treated as discrete.  May contain entries that do not appear in
#'     colnames(x)
#' @param discrete_var_range_fns a list with one entry for each element
#'     of discrete_vars.  Each entry is a named list of length 2; the element
#'     named "a" is a character string with the name of a function that returns
#'     a(x) for any real x, and the element named "b" is a character string with
#'     the name of a function that returns b(x) for any real x.
#' @param lower Vector of lower truncation points
#' @param upper Vector of upper truncation points
#' @param log logical; if TRUE, return the log of the kernel function value
#' @param ... mop up extra arguments
#' 
#' @return the value of the kernel function given by the pdtmvn distribution
#'     at x.
log_pdtmvn_mode_centered_kernel <- function(x,
		center,
		bw,
		bw_continuous,
		conditional_bw_discrete,
		conditional_center_discrete_offset_multiplier,
		continuous_vars,
		discrete_vars,
		continuous_var_col_inds,
		discrete_var_col_inds,
		discrete_var_range_fns,
		lower,
		upper,
        x_names,
		log,
        validate_in_support = TRUE,
        validate_level = 1L,
        ...) {
    if(missing(bw_continuous)) {
        bw_continuous <- NULL
    }
    if(missing(conditional_bw_discrete)) {
        conditional_bw_discrete <- NULL
    }
    if(missing(conditional_center_discrete_offset_multiplier)) {
        conditional_center_discrete_offset_multiplier <- NULL
    }
    
    ## center parameter of pdtmvn_kernel is mean of log
    ## mode of resulting log-normal distribution is
    ## mode = exp(mu - bw %*% 1) (where 1 is a column vector of 1s)
    ## therefore mu = log(mode) + bw %*% 1
    ## we have to jump through some hoops to ensure that the
    ## variables come in the right order
    reduced_x_names <- names(x)
    inds_x_vars_in_orig_vars <- which(x_names %in% reduced_x_names)
    x_names_for_call <- x_names[inds_x_vars_in_orig_vars]
    
    mean_offset <- apply(bw, 1, sum)[x_names %in% colnames(center)]
    unadjusted_result <- pdtmvn_kernel(x = log(x)[, x_names_for_call, drop = FALSE],
        center = sweep(log(center)[, x_names_for_call, drop = FALSE], 2, mean_offset, `+`),
        bw = bw,
        bw_continuous = bw_continuous,
        conditional_bw_discrete = conditional_bw_discrete,
        conditional_center_discrete_offset_multiplier = conditional_center_discrete_offset_multiplier,
        continuous_vars = continuous_vars,
        discrete_vars = discrete_vars,
        continuous_var_col_inds = continuous_var_col_inds,
        discrete_var_col_inds = discrete_var_col_inds,
        discrete_var_range_fns = discrete_var_range_fns,
        lower = lower,
        upper = upper,
        x_names = x_names,
        log = log,
        validate_in_support = validate_in_support,
        validate_level = validate_level
    )
    
    continuous_vars_in_x <- continuous_vars[continuous_vars %in% colnames(x)]
    if(length(continuous_vars_in_x) > 0) {
        log_adjustment_factor <- -1 * apply(log(x[, continuous_vars_in_x, drop = FALSE]), 1, sum)
    } else {
        log_adjustment_factor <- rep(0, nrow(x))
    }
    
    if(log) {
        return(unadjusted_result + log_adjustment_factor)
    } else {
        return(unadjusted_result * exp(log_adjustment_factor))
    }
}

#' Simulate from the kernel function given by the pdtmvn distribution.
#' 
#' @param n number of simulations to generate
#' @param a matrix of values at which to evaluate the kernel function, with
#'     column names specified.  Each row is an observation, each column is an
#'     observed variable.
#' @param center a real vector, center point for the kernel function
#' @param bw bandwidth matrix
#' @param bw_continuous the portion of bw corresponding to continuous variables
#' @param conditional_bw_discrete the Schur complement of the portion of bw
#'     corresponding to discrete variables
#' @param conditional_center_discrete_offset_multiplier Sigma_dc Sigma_c^{-1}.
#'     This is used in computing the mean of the underlying multivariate normal
#'     distribution for the discrete variables conditioning on the continuous
#'     variables.
#' @param continuous_vars character vector with names of variables that are to
#'     be treated as continuous.  May contain entries that do not appear in
#'     colnames(x).
#' @param discrete_vars character vector with names of variables that are to be
#'     treated as discrete.  May contain entries that do not appear in
#'     colnames(x)
#' @param discrete_var_range_fns a list with one entry for each element
#'     of discrete_vars.  Each entry is a named list of length 2; the element
#'     named "a" is a character string with the name of a function that returns
#'     a(x) for any real x, and the element named "b" is a character string with
#'     the name of a function that returns b(x) for any real x.
#' @param lower Vector of lower truncation points
#' @param upper Vector of upper truncation points
#' @param log logical; if TRUE, return the log of the kernel function value
#' @param ... mop up extra arguments
#' 
#' @return the value of the kernel function given by the pdtmvn distribution
#'     at x.
rlog_pdtmvn_mode_centered_kernel <- function(n,
    conditioning_obs,
    center,
    bw,
    bw_continuous,
    conditional_bw_discrete,
    conditional_center_discrete_offset_multiplier,
    continuous_vars,
    discrete_vars,
    continuous_var_col_inds,
    discrete_var_col_inds,
    discrete_var_range_fns,
    lower,
    upper,
    x_names,
    ...) {
    if(missing(conditioning_obs) || is.null(conditioning_obs)) {
        log_conditioning_obs <- NULL
    } else {
        log_conditioning_obs <- log(conditioning_obs)
    }
    if(missing(bw_continuous)) {
        bw_continuous <- NULL
    }
    if(missing(conditional_bw_discrete)) {
        conditional_bw_discrete <- NULL
    }
    if(missing(conditional_center_discrete_offset_multiplier)) {
        conditional_center_discrete_offset_multiplier <- NULL
    }
    
    ## center parameter of pdtmvn_kernel is mean of log
    ## mode of resulting log-normal distribution is
    ## mode = exp(mu - bw %*% 1) (where 1 is a column vector of 1s)
    ## therefore mu = log(mode) + bw %*% 1
    reduced_x_names <- names(center)
    inds_x_vars_in_orig_vars <- which(x_names %in% reduced_x_names)
    x_names_for_call <- x_names[inds_x_vars_in_orig_vars]
    
    mean_offset <- apply(bw, 1, sum)[x_names %in% colnames(center)]
        
    return(exp(rpdtmvn_kernel(n = n,
                conditioning_obs = log_conditioning_obs,
                center = sweep(log(center)[, x_names_for_call, drop = FALSE], 2, mean_offset, `+`),
                bw = bw,
                bw_continuous = bw_continuous,
                conditional_bw_discrete = conditional_bw_discrete,
                conditional_center_discrete_offset_multiplier = conditional_center_discrete_offset_multiplier,
                continuous_vars = continuous_vars,
                discrete_vars = discrete_vars,
                continuous_var_col_inds = continuous_var_col_inds,
                discrete_var_col_inds = discrete_var_col_inds,
                discrete_var_range_fns = discrete_var_range_fns,
                lower = lower,
                upper = upper,
                x_names = x_names)[, reduced_x_names, drop = FALSE]))
}


#' A function to vectorize the parameters of the pdtmvn_kernel and convert
#' to estimation scale.
#' 
#' @param theta_list parameters for the pdtmvn_kernel in list format
#' @param parameterization character; currently, only supported value is
#'     "bw-diagonalized-est-eigenvalues"
#' @param ... mop up other arguments
#'
#' @return vector containing parameters that are estimated on a scale
#'     suitable for numerical optimization
vectorize_params_log_pdtmvn_mode_centered_kernel <- function(theta_list, ...) {
    return(vectorize_params_pdtmvn_kernel(theta_list = theta_list))
}

#' A function to unvectorize the parameters of the pdtmvn_kernel and convert
#' from estimation scale to scale actually used.
#' 
#' @param theta_vector a numeric vector of parameters that are being optimized,
#'     on a scale suitable for use in optim.
#' @param parameterization character; currently, only supported value is
#'     "bw-diagonalized-est-eigenvalues"
#' @param bw_evecs square matrix with eigenvectors of the bandwidth matrix in
#'     columns
#' @param continuous_vars character vector with variables names that are to be
#'     treated as continuous
#' @param discrete_vars character vector with variables names that are to be
#'     treated as discrete
#' @param kcde_control list of control parameters to kcde
#' 
#' @return list of parameters to pdtmvn_kernel
update_theta_from_vectorized_theta_est_log_pdtmvn_mode_centered_kernel <- function(theta_est_vector, theta) {
    return(update_theta_from_vectorized_theta_est_pdtmvn_kernel(theta_est_vector = theta_est_vector, theta = theta))
}

#' Get initial parameter values for the pdtmvn_kernel
#' 
#' @param parameterization Character vector of length 1; the parameterization to
#'     use.  Currently, the only supported option is
#'     "bw-diagonalized-est-eigenvalues"
#' @param x a matrix of values used to initialize the parameter values for the
#'     kernel function.  Each row is an observation, each column is an
#'     observed variable.
#' @param continuous_vars character vector with variables names that are to be
#'     treated as continuous
#' @param discrete_vars character vector with variables names that are to be
#'     treated as discrete
#' @param kcde_control list of control parameters to kcde
#' @param ... used to absorb other arguments in the function call
#' 
#' @return list with initial values of parameters to the pdtmvn_kernel
initialize_params_log_pdtmvn_mode_centered_kernel <- function(prev_theta,
	x,
    total_num_vars,
    sample_size,
	...) {
	
    return(initialize_params_pdtmvn_kernel(prev_theta = prev_theta,
    	x = log(x),
        total_num_vars = total_num_vars,
        sample_size = sample_size))
}


#' Get lower and upper bounds for the theta parameters being estimated
#' in the pdtmvn_kernel
#' 
#' @param theta_list parameters to pdtmvn kernel in list format
#' @param ... mop up arguments
#' 
#' @return list with two components: lower and upper, numeric vectors
get_theta_optim_bounds_log_pdtmvn_mode_centered_kernel <- function(theta_list, ...) {
    return(get_theta_optim_bounds_pdtmvn_kernel(theta_list))
}
