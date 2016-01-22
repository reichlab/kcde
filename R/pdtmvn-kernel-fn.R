## pdtmvn kernel function
## 
## get_col_inds_continuous_dicrete_vars_used
## pdtmvn_kernel
## rpdtmvn_kernel
## compute_pdtmvn_kernel_bw_params_from_bw_eigen
## vectorize_params_pdtmvn_kernel
## update_theta_from_vectorized_theta_est_pdtmvn_kernel
## initialize_params_pdtmvn_kernel

#' Create two integer vectors with indices of columns in x corresponding to
#' continuous variables and discrete variables.  These integer vectors are
#' suitable for passing to functions from the pdtmvn function.
#' 
#' @param x_colnames character vector with column names of x
#' @param continuous_vars character vector with names of columns in x that
#'     should be treated as continuous variables
#' @param discrete_vars character vector with names of columns in x that
#'     should be treated as discrete variables
#' 
#' @return list with two components: continuous_vars is an integer vector of
#'     columns in x corresponding to continuous variables and discrete_vars is
#'     an integer vector of columns in x corresponding to discrete variables
get_col_inds_continuous_discrete_vars_used <- function(x_colnames,
		continuous_vars,
		discrete_vars) {
	if(is.null(x_colnames)) {
		stop("x must have column names")
	}
	col_inds_continuous_vars <- seq_along(x_colnames)[x_colnames %in% continuous_vars]
	col_inds_discrete_vars <- seq_along(x_colnames)[x_colnames %in% discrete_vars]
	
	return(list(
		continuous_vars = col_inds_continuous_vars,
		discrete_vars = col_inds_discrete_vars
	))
}

#' Evaluate the kernel function given by the pdtmvn distribution.
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
pdtmvn_kernel <- function(x,
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
        ...) {
    if(length(dim(x)) > 0) {
        x_arg_names <- colnames(x)
        x <- as.vector(as.matrix(x))
        names(x) <- x_arg_names
    }
    
    if(ncol(center) != ncol(bw)) {
        ## evaluating a marginal kernel for a subset of variables
        if(any(c(lower != -Inf, upper != Inf)) && ncol(center) != ncol(bw)) {
            stop("Lower and upper bounds must be -Inf and Inf respectively if evaluating a marginal kernel.")
        }
        
        
        inds_x_vars_in_orig_vars <- which(x_names %in% names(x))
        
        return(pdtmvn::dpdtmvn(x = center[, x_names[inds_x_vars_in_orig_vars]],
                mean = x[x_names[inds_x_vars_in_orig_vars]],
                sigma = bw[inds_x_vars_in_orig_vars, inds_x_vars_in_orig_vars, drop = FALSE],
#                sigma_continuous = bw_continuous,
#                conditional_sigma_discrete = conditional_bw_discrete,
#                conditional_mean_discrete_offset_multiplier = 
#                    conditional_center_discrete_offset_multiplier,
                lower = lower[names(x)],
                upper = upper[names(x)],
#                continuous_vars = continuous_var_col_inds,
#                discrete_vars = discrete_var_col_inds,
                discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
                log = log,
                validate_level = 1L))
    } else {
        return(pdtmvn::dpdtmvn(x = center[, x_names],
                mean = x[x_names],
                sigma = bw,
                sigma_continuous = bw_continuous,
                conditional_sigma_discrete = conditional_bw_discrete,
                conditional_mean_discrete_offset_multiplier = 
                    conditional_center_discrete_offset_multiplier,
                lower = lower[names(x)],
                upper = upper[names(x)],
                continuous_vars = continuous_var_col_inds,
                discrete_vars = discrete_var_col_inds,
                discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
                log = log,
                validate_level = 1L))
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
rpdtmvn_kernel <- function(n,
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
    if(length(dim(center)) > 0) {
        center_names <- colnames(center)
        center <- as.vector(as.matrix(center))
        names(center) <- center_names
    }
    return(pdtmvn::rpdtmvn(n = n,
            x_fixed = conditioning_obs,
            mean = t(as.matrix(center))[, x_names],
            sigma = bw,
            lower = lower[x_names],
            upper = upper[x_names],
            continuous_vars = continuous_var_col_inds,
            discrete_vars = discrete_var_col_inds,
            discrete_var_range_fns = discrete_var_range_fns,
            validate_level = 1))
}


#' Compute the parameters bw, bw_continuous, conditional_bw_discrete, and
#' conditional_center_discrete_offset_multiplier for the pdtmvn_kernel function
#' from the eigen-decomposition of the bandwidth matrix.
#' 
#' @param bw_evecs matrix whose columns are the eigenvectors of the bandwidth
#'     matrix.
#' @param bw_evals vector of eigenvalues of the bandwidth matrix.
#' @param continuous_vars Vector containing column indices for continuous
#'     variables.
#' @param discrete_vars Vector containing column indices for discrete variables.
#' 
#' @return Named list with four components: bw, bw_continuous,
#'     conditional_bw_discrete, and conditional_center_discrete_offset_multiplier
compute_pdtmvn_kernel_bw_params_from_bw_eigen <- function(bw_evecs,
    bw_evals,
    continuous_var_col_inds,
    discrete_var_col_inds) {
    bw <- bw_evecs %*% diag(bw_evals, nrow = length(bw_evals)) %*% t(bw_evecs)
    
    bw_params <- c(
		list(bw = bw,
			bw_evecs = bw_evecs,
			bw_evals = bw_evals,
			log_bw_evals = log(bw_evals)),
        pdtmvn::compute_sigma_subcomponents(sigma = bw,
            continuous_vars = continuous_var_col_inds,
            discrete_vars = discrete_var_col_inds,
            validate_level = 0)
    )
    
    names(bw_params)[names(bw_params) == "sigma_continuous"] <- "bw_continuous"
    names(bw_params)[names(bw_params) == "conditional_sigma_discrete"] <- 
        "conditional_bw_discrete"
    names(bw_params)[names(bw_params) == "conditional_mean_discrete_offset_multiplier"] <- 
        "conditional_center_discrete_offset_multiplier"
    
    return(bw_params)
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
vectorize_params_pdtmvn_kernel <- function(theta_list, ...) {
	if(identical(theta_list$parameterization, "bw-diagonalized-est-eigenvalues")) {
		return(theta_list$log_bw_evals)
	} else if(identical(theta_list$parameterization, "bw-chol-decomp")) {
        return(theta_list$bw_chol_decomp_vec)
    } else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
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
update_theta_from_vectorized_theta_est_pdtmvn_kernel <- function(theta_est_vector, theta) {
	if(identical(theta$parameterization, "bw-diagonalized-est-eigenvalues")) {
		num_bw_evals <- ncol(theta$bw_evecs)
		temp <- compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = theta$bw_evecs,
			bw_evals = exp(theta_est_vector[seq_len(num_bw_evals)]),
			theta$continuous_var_col_inds,
			theta$discrete_var_col_inds)
		
		theta[names(temp)] <- temp

		return(list(
			theta = theta,
			num_theta_vals_used = num_bw_evals
		))
    } else if(identical(theta_list$parameterization, "bw-chol-decomp")) {
        num_bw_chol_decomp_vec_entries <- ncol(theta$bw_evecs) * (1 + ncol(theta$bw_evecs)) / 2
        temp <- compute_pdtmvn_kernel_bw_params_from_bw_chol_decomp(bw_chol_decomp = ,
            theta$continuous_var_col_inds,
            theta$discrete_var_col_inds)
        
        theta[names(temp)] <- temp
        
        return(list(
            theta = theta,
            num_theta_vals_used = num_bw_evals
        ))
    } else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
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
initialize_params_pdtmvn_kernel <- function(prev_theta,
	x,
	...) {
	
    new_theta <- prev_theta
    
	if(identical(new_theta$parameterization, "bw-diagonalized-est-eigenvalues")) {
		continuous_discrete_var_col_inds <-
			get_col_inds_continuous_discrete_vars_used(colnames(x),
                new_theta$continuous_vars,
                new_theta$discrete_vars)
		
        if(ncol(x) > 1) {
            new_theta$x_names <- colnames(x)
            
            require(robust)
            sample_cov_hat <- robust::covRob(x)$cov
        } else {
            new_theta$x_names <- names(x)
            
            sample_cov_hat <- matrix(var(x))
        }
		sample_cov_eigen <- eigen(sample_cov_hat)
		
        bw_params <- compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = sample_cov_eigen$vectors,
		    bw_evals = sample_cov_eigen$values,
			continuous_discrete_var_col_inds$continuous_vars,
			continuous_discrete_var_col_inds$discrete_vars)
		new_theta[names(bw_params)] <- bw_params
        
		new_theta$continuous_var_col_inds = continuous_discrete_var_col_inds$continuous_vars
        new_theta$discrete_var_col_inds = continuous_discrete_var_col_inds$discrete_vars
	} else {
		stop("Invalid parameterization for pdtmvn kernel function")
	}
    
    return(new_theta)
}
