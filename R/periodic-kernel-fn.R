## periodic kernel function
## 
## periodic_kernel

#' Evaluate the periodic kernel function
#' 
#' @param x a vector of values at which to evaluate the kernel function
#' @param center a real number, center point for the kernel function
#' @param period kernel period
#' @param bw kernel bandwidth
#' 
#' @return vector of the kernel function value at each point in x
periodic_kernel <- function(x, center, period, bw, log, ...) {
    result <- -0.5 * (sin(period * (as.vector(as.matrix(x)) - as.vector(as.matrix(center)))) / bw)^2
    
    if(log) {
        return(result)
    } else {
        return(exp(result))
    }
}

#' Get lower and upper bounds for the log_bw parameter in the periodic_kernel
#' 
#' @param ... mop up arguments
#' 
#' @return list with two components: lower and upper, numeric vectors
get_theta_optim_bounds_periodic_kernel <- function(...) {
    return(list(
        lower = -50,
        upper = Inf
    ))
}

#' A function to vectorize the parameters of the pdtmvn_kernel and convert
#' to estimation scale.
#' @param theta_list parameters for the pdtmvn_kernel in list format
#' @param parameterization character; currently, only supported value is
#'     "bw-diagonalized-est-eigenvalues"
#' @param ... mop up other arguments
#'
#' @return vector containing parameters that are estimated on a scale
#'     suitable for numerical optimization
vectorize_params_periodic_kernel <- function(theta_list, parameterization, ...) {
    return(theta_list$log_bw)
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
update_theta_from_vectorized_theta_est_periodic_kernel <- function(theta_est_vector, theta) {
	theta$log_bw <- theta_est_vector[1]
    theta$bw <- exp(theta$log_bw)
    
    return(list(
		theta = theta,
		num_theta_vals_used = 1L
	))
}


#' Get initial parameter values for the periodic_kernel
#' 
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
#' @return list with initial values of parameters to the periodic_kernel
initialize_params_periodic_kernel <- function(x,
    prev_theta,
    total_num_vars,
    sample_size,
    ...) {
    new_theta <- prev_theta
    new_theta$bw <- sample_size^{-2/(total_num_vars + 4)} * sqrt(0.5)
    new_theta$log_bw <- log(new_theta$bw)
    
    return(new_theta)
}
