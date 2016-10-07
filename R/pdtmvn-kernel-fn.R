## pdtmvn kernel function
##
## get_col_inds_continuous_discrete_vars_used (exported)
## pdtmvn_kernel (exported)
## rpdtmvn_kernel (exported)
## compute_pdtmvn_kernel_bw_params_from_bw_eigen (exported)
## compute_pdtmvn_kernel_bw_params_from_bw_chol_decomp (exported)
## get_theta_optim_bounds_pdtmvn_kernel (exported)
## vectorize_params_pdtmvn_kernel (exported)
## update_theta_from_vectorized_theta_est_pdtmvn_kernel (exported)
## initialize_params_pdtmvn_kernel (exported)

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
#' @export
get_col_inds_continuous_discrete_vars_used <- function(x_colnames,
                                                       continuous_vars,
                                                       discrete_vars) {
    if (is.null(x_colnames)) {
        stop("x must have column names")
    }
    col_inds_continuous_vars <-
        seq_along(x_colnames)[x_colnames %in% continuous_vars]
    col_inds_discrete_vars <-
        seq_along(x_colnames)[x_colnames %in% discrete_vars]

    return(
        list(continuous_vars = col_inds_continuous_vars,
             discrete_vars = col_inds_discrete_vars)
    )
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
#' @export
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
                          validate_in_support = TRUE,
                          validate_level = 1L,
                          ...) {
    ## if(length(dim(x)) > 0) {
    ##      x_arg_names <- colnames(x)
    ##      x <- as.vector(as.matrix(x))
    ##      names(x) <- x_arg_names
    ## }
    if (length(dim(center)) == 0) {
        center_arg_names <- names(center)
        dim(center) <- c(1, length(center))
        colnames(center) <- center_arg_names
    } else if (!is.matrix(center)) {
        center <- as.matrix(center)
    }

    if (ncol(x) != ncol(bw)) {
        ## evaluating a marginal kernel for a subset of variables
        if (any(c(lower != -Inf, upper != Inf)) &&
            ncol(center) != ncol(bw)) {
            stop(
                "Lower and upper bounds must be -Inf and Inf respectively if evaluating a marginal kernel."
            )
        }

        reduced_x_names <- colnames(x)
        inds_x_vars_in_orig_vars <-
            which(x_names %in% reduced_x_names)
        continuous_x_names <- x_names[continuous_var_col_inds]
        reduced_continuous_var_col_inds <-
            which(reduced_x_names %in% continuous_x_names)
        discrete_x_names <- x_names[discrete_var_col_inds]
        reduced_discrete_var_col_inds <-
            which(reduced_x_names %in% discrete_x_names)

        x_names_for_call <- x_names[inds_x_vars_in_orig_vars]
        continuous_var_col_inds_for_call <-
            reduced_continuous_var_col_inds
        discrete_var_col_inds_for_call <-
            reduced_discrete_var_col_inds

        ## return(sapply(seq_len(nrow(center)), function(center_row_ind) {
        ##     pdtmvn::dpdtmvn(x = x[, x_names_for_call, drop = FALSE],
        ##                     mean = center[center_row_ind, x_names[inds_x_vars_in_orig_vars]],
        ##                     sigma = bw[inds_x_vars_in_orig_vars, inds_x_vars_in_orig_vars, drop = FALSE],
        ##                     sigma_continuous = bw_continuous,
        ##                     conditional_sigma_discrete = conditional_bw_discrete,
        ##                     conditional_mean_discrete_offset_multiplier =
        ##                         conditional_center_discrete_offset_multiplier,
        ##                     lower = lower[colnames(x)],
        ##                     upper = upper[colnames(x)],
        ##                     continuous_vars = reduced_continuous_var_col_inds,
        ##                     discrete_vars = reduced_discrete_var_col_inds,
        ##                     discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
        ##                     log = log,
        ##                     validate_in_support = validate_in_support,
        ##                     validate_level = validate_level)
        ## }))
    } else {
        x_names_for_call <- x_names
        inds_x_vars_in_orig_vars <- seq_len(nrow(bw))
        continuous_var_col_inds_for_call <- continuous_var_col_inds
        discrete_var_col_inds_for_call <- discrete_var_col_inds
        ## return(sapply(seq_len(nrow(center)), function(center_row_ind) {
        ##     pdtmvn::dpdtmvn(x = x[, x_names, drop = FALSE],
        ##                     mean = center[center_row_ind, x_names],
        ##                     sigma = bw,
        ##                     sigma_continuous = bw_continuous,
        ##                     conditional_sigma_discrete = conditional_bw_discrete,
        ##                     conditional_mean_discrete_offset_multiplier =
        ##                         conditional_center_discrete_offset_multiplier,
        ##                     lower = lower[colnames(x)],
        ##                     upper = upper[colnames(x)],
        ##                     continuous_vars = continuous_var_col_inds,
        ##                     discrete_vars = discrete_var_col_inds,
        ##                     discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
        ##                     log = log,
        ##                     validate_in_support = validate_in_support,
        ##                     validate_level = validate_level)
        ## }))
    }

    if (nrow(center) > 1) {
        return(sapply(seq_len(nrow(center)), function(center_row_ind) {
            pdtmvn::dpdtmvn(
                x = x[, x_names_for_call, drop = FALSE],
                mean = center[center_row_ind, x_names[inds_x_vars_in_orig_vars]],
                sigma = bw[inds_x_vars_in_orig_vars, inds_x_vars_in_orig_vars, drop = FALSE],
                ## sigma_continuous = bw_continuous,
                ## conditional_sigma_discrete = conditional_bw_discrete,
                ## conditional_mean_discrete_offset_multiplier =
                ##     conditional_center_discrete_offset_multiplier,
                lower = lower[colnames(x)],
                upper = upper[colnames(x)],
                continuous_vars = continuous_var_col_inds_for_call,
                discrete_vars = discrete_var_col_inds_for_call,
                discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
                log = log,
                validate_in_support = validate_in_support,
                validate_level = validate_level
            )
        }))
    } else {
        return(
            pdtmvn::dpdtmvn(
                x = x[, x_names_for_call, drop = FALSE],
                mean = center[1, ],
                sigma = bw[inds_x_vars_in_orig_vars, inds_x_vars_in_orig_vars, drop = FALSE],
                ## sigma_continuous = bw_continuous,
                ## conditional_sigma_discrete = conditional_bw_discrete,
                ## conditional_mean_discrete_offset_multiplier =
                ##     conditional_center_discrete_offset_multiplier,
                lower = lower[colnames(x)],
                upper = upper[colnames(x)],
                continuous_vars = continuous_var_col_inds_for_call,
                discrete_vars = discrete_var_col_inds_for_call,
                discrete_var_range_fns = discrete_var_range_fns[discrete_vars],
                log = log,
                validate_in_support = validate_in_support,
                validate_level = validate_level
            )
        )
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
#' @export
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
    ## if(length(dim(center)) > 0) {
    ##     center_names <- colnames(center)
    ##     center <- as.vector(as.matrix(center))
    ##     names(center) <- center_names
    ## }
    if (length(dim(center)) == 0) {
        center_arg_names <- names(center)
        dim(center) <- c(1, length(center))
        colnames(center) <- center_arg_names
    } else if (!identical(class(center), "matrix")) {
        center <- as.matrix(center)
    }

    if (missing(conditioning_obs)) {
        conditioning_obs <- NULL
    }
    return(
        pdtmvn::rpdtmvn(
            n = n,
            x_fixed = conditioning_obs,
            mean = center[, x_names, drop = FALSE],
            sigma = bw,
            lower = lower[x_names],
            upper = upper[x_names],
            continuous_vars = continuous_var_col_inds,
            discrete_vars = discrete_var_col_inds,
            discrete_var_range_fns = discrete_var_range_fns,
            validate_level = 1
        )
    )
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
#' @export
compute_pdtmvn_kernel_bw_params_from_bw_eigen <-
    function(multilogit_bw_evecs,
             bw_evecs,
             log_bw_evals,
             bw_evals,
             continuous_var_col_inds,
             discrete_var_col_inds) {
        bw <-
            bw_evecs %*% diag(bw_evals, nrow = length(bw_evals)) %*% t(bw_evecs)

        bw_params <- c(
            list(
                bw = bw,
                multilogit_bw_evecs = multilogit_bw_evecs,
                bw_evecs = bw_evecs,
                log_bw_evals = log_bw_evals,
                bw_evals = bw_evals
            ),
            pdtmvn::compute_sigma_subcomponents(
                sigma = bw,
                continuous_vars = continuous_var_col_inds,
                discrete_vars = discrete_var_col_inds,
                validate_level = 0
            )
        )

        names(bw_params)[names(bw_params) == "sigma_continuous"] <-
            "bw_continuous"
        names(bw_params)[names(bw_params) == "conditional_sigma_discrete"] <-
            "conditional_bw_discrete"
        names(bw_params)[names(bw_params) == "conditional_mean_discrete_offset_multiplier"] <-
            "conditional_center_discrete_offset_multiplier"

        return(bw_params)
    }

#' Compute the parameters bw, bw_continuous, conditional_bw_discrete, and
#' conditional_center_discrete_offset_multiplier for the pdtmvn_kernel function
#' from the Cholesky decomposition of the bandwidth matrix.
#'
#' @param bw_chol_decomp lower triangular matrix giving the Cholesky
#'     decomposition of the bandwidth matrix.
#' @param continuous_vars Vector containing column indices for continuous
#'     variables.
#' @param discrete_vars Vector containing column indices for discrete variables.
#'
#' @return Named list with four components: bw, bw_continuous,
#'     conditional_bw_discrete, and conditional_center_discrete_offset_multiplier
#' @export
compute_pdtmvn_kernel_bw_params_from_bw_chol_decomp <-
    function(bw_chol_decomp,
             bw_chol_decomp_vec,
             continuous_var_col_inds,
             discrete_var_col_inds) {
        bw <- bw_chol_decomp %*% t(bw_chol_decomp)

        bw_params <- c(
            list(
                bw = bw,
                bw_chol_decomp = bw_chol_decomp,
                bw_chol_decomp_vec = bw_chol_decomp_vec
            ),
            pdtmvn::compute_sigma_subcomponents(
                sigma = bw,
                continuous_vars = continuous_var_col_inds,
                discrete_vars = discrete_var_col_inds,
                validate_level = 0
            )
        )

        names(bw_params)[names(bw_params) == "sigma_continuous"] <-
            "bw_continuous"
        names(bw_params)[names(bw_params) == "conditional_sigma_discrete"] <-
            "conditional_bw_discrete"
        names(bw_params)[names(bw_params) == "conditional_mean_discrete_offset_multiplier"] <-
            "conditional_center_discrete_offset_multiplier"

        return(bw_params)
    }

#' Get lower and upper bounds for the theta parameters being estimated
#' in the pdtmvn_kernel
#'
#' @param theta_list parameters to pdtmvn kernel in list format
#' @param ... mop up arguments
#'
#' @return list with two components: lower and upper, numeric vectors
#' @export
get_theta_optim_bounds_pdtmvn_kernel <- function(theta_list, ...) {
    if (identical(theta_list$parameterization,
                  "bw-diagonalized-est-eigenvalues")) {
        return(list(lower = -50,
                    upper = Inf))
    } else if (identical(theta_list$parameterization, "bw-diagonalized")) {
        return(list(lower = -50,
                    upper = Inf))
    } else if (identical(theta_list$parameterization, "bw-chol-decomp")) {
        lower_triang_inds <- as.matrix(expand.grid(seq_len(
            ncol(theta_list$bw_chol_decomp)
        ),
        seq_len(
            ncol(theta_list$bw_chol_decomp)
        )))
        lower_triang_inds <-
            lower_triang_inds[lower_triang_inds[, 1] >= lower_triang_inds[, 2], , drop = FALSE]
        lower_bounds_matrix <-
            matrix(-50,
                   nrow = nrow(theta_list$bw),
                   ncol = ncol(theta_list$bw))
        diag(lower_bounds_matrix) <- 10 ^ {
            -6
        }
        upper_bounds_matrix <-
            matrix(Inf,
                   nrow = nrow(theta_list$bw),
                   ncol = ncol(theta_list$bw))
        return(list(lower = lower_bounds_matrix[lower_triang_inds],
                    upper = upper_bounds_matrix[lower_triang_inds]))
    } else {
        stop("Invalid parameterization for pdtmvn kernel function")
    }
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
#' @export
vectorize_params_pdtmvn_kernel <- function(theta_list, ...) {
    if (identical(theta_list$parameterization,
                  "bw-diagonalized-est-eigenvalues")) {
        return(theta_list$log_bw_evals)
    } else if (identical(theta_list$parameterization, "bw-diagonalized")) {
        return(c(
            theta_list$log_bw_evals,
            theta_list$multilogit_bw_evecs
        ))
    } else if (identical(theta_list$parameterization, "bw-chol-decomp")) {
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
#' @export
update_theta_from_vectorized_theta_est_pdtmvn_kernel <-
    function(theta_est_vector, theta) {
        if (identical(theta$parameterization,
                      "bw-diagonalized-est-eigenvalues")) {
            num_bw_evals <- ncol(theta$bw_evecs)
            temp <- compute_pdtmvn_kernel_bw_params_from_bw_eigen(
                multilogit_bw_evecs = theta$multilogit_bw_evecs,
                bw_evecs = theta$bw_evecs,
                log_bw_evals = theta_est_vector[seq_len(num_bw_evals)],
                bw_evals = exp(theta_est_vector[seq_len(num_bw_evals)]),
                theta$continuous_var_col_inds,
                theta$discrete_var_col_inds
            )

            theta[names(temp)] <- temp

            return(list(theta = theta,
                        num_theta_vals_used = num_bw_evals))
        } else if (identical(theta$parameterization, "bw-diagonalized")) {
            num_bw_evals <- ncol(theta$bw_evecs)
            num_bw_evec_params <-
                ncol(theta$bw_evecs) * (ncol(theta$bw_evecs) - 1)
            multilogit_bw_evecs <-
                theta_est_vector[num_bw_evals + seq_len(num_bw_evec_params)]
            updated_evecs <-
                compute_bw_evecs_from_multilogit_bw_evecs(multilogit_bw_evecs,
                                                          ncol(theta$bw_evecs))
            temp <- compute_pdtmvn_kernel_bw_params_from_bw_eigen(
                multilogit_bw_evecs = multilogit_bw_evecs,
                bw_evecs = theta$bw_evecs,
                log_bw_evals = theta_est_vector[seq_len(num_bw_evals)],
                bw_evals = exp(theta_est_vector[seq_len(num_bw_evals)]),
                theta$continuous_var_col_inds,
                theta$discrete_var_col_inds
            )

            theta[names(temp)] <- temp

            return(list(
                theta = theta,
                num_theta_vals_used = num_bw_evals + num_bw_evec_params
            ))
        } else if (identical(theta$parameterization, "bw-chol-decomp")) {
            ## Obtain vector of parameters for Cholesky decomposition
            num_bw_chol_decomp_vec_entries <-
                ncol(theta$bw_chol_decomp) * (1 + ncol(theta$bw_chol_decomp)) / 2
            bw_chol_decomp_vec <-
                theta_est_vector[seq_len(num_bw_chol_decomp_vec_entries)]

            ## Update Cholesky decomposition matrix
            update_inds <- as.matrix(expand.grid(seq_len(ncol(
                theta$bw_chol_decomp
            )),
            seq_len(ncol(
                theta$bw_chol_decomp
            ))))
            update_inds <-
                update_inds[update_inds[, 1] >= update_inds[, 2], , drop = FALSE]
            theta$bw_chol_decomp[update_inds] <- bw_chol_decomp_vec

            ## Enforce lower bound for entries on diagonal:
            ## > 0 by an amount that's large enough to prevent numerical instability
            zero_threshold <- 10 ^ -7
            diag(theta$bw_chol_decomp)[diag(theta$bw_chol_decomp) < zero_threshold] <-
                zero_threshold

            ## Compute other bandwidth parameters used in calls to pdtmvn kernel
            temp <-
                compute_pdtmvn_kernel_bw_params_from_bw_chol_decomp(
                    bw_chol_decomp_vec = bw_chol_decomp_vec,
                    bw_chol_decomp = theta$bw_chol_decomp,
                    theta$continuous_var_col_inds,
                    theta$discrete_var_col_inds
                )

            theta[names(temp)] <- temp

            ## Return
            return(list(theta = theta,
                        num_theta_vals_used = num_bw_chol_decomp_vec_entries))
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
#' @export
initialize_params_pdtmvn_kernel <- function(prev_theta,
                                            x,
                                            total_num_vars,
                                            sample_size,
                                            ...) {
    new_theta <- prev_theta

    if (new_theta$parameterization %in% c("bw-diagonalized", "bw-diagonalized-est-eigenvalues")) {
        continuous_discrete_var_col_inds <-
            get_col_inds_continuous_discrete_vars_used(colnames(x),
                                                       new_theta$continuous_vars,
                                                       new_theta$discrete_vars)

        if (ncol(x) > 1) {
            new_theta$x_names <- colnames(x)

            require(robust)
            sample_cov_hat <-
                tryCatch(
                    robust::covRob(x)$cov,
                    error = function(e) {
                        cov(x)
                    }
                )
        } else {
            new_theta$x_names <- names(x)

            sample_cov_hat <- matrix(var(x))
        }
        init_bw <-
            sample_size ^ {
                -2 / (total_num_vars + 4)
            } * sample_cov_hat
        init_bw_eigen <- eigen(init_bw)

        bw_params <- compute_pdtmvn_kernel_bw_params_from_bw_eigen(
            multilogit_bw_evecs = compute_multilogit_bw_evecs_from_bw_evecs(init_bw_eigen$vectors),
            bw_evecs = init_bw_eigen$vectors,
            log_bw_evals = log(init_bw_eigen$values),
            bw_evals = init_bw_eigen$values,
            continuous_discrete_var_col_inds$continuous_vars,
            continuous_discrete_var_col_inds$discrete_vars
        )
        new_theta[names(bw_params)] <- bw_params

        new_theta$continuous_var_col_inds = continuous_discrete_var_col_inds$continuous_vars
        new_theta$discrete_var_col_inds = continuous_discrete_var_col_inds$discrete_vars
    } else if (identical(new_theta$parameterization, "bw-chol-decomp")) {
        continuous_discrete_var_col_inds <-
            get_col_inds_continuous_discrete_vars_used(colnames(x),
                                                       new_theta$continuous_vars,
                                                       new_theta$discrete_vars)

        if (ncol(x) > 1) {
            new_theta$x_names <- colnames(x)

            require(robust)
            sample_cov_hat <-
                tryCatch(
                    robust::covRob(x)$cov,
                    error = function(e) {
                        cov(x)
                    }
                )
        } else {
            new_theta$x_names <- names(x)

            sample_cov_hat <- matrix(var(x))
        }
        ## Get Cholesky decomposition.  R returns upper triangular portion,
        ## but our parameterization works with lower triangular portion, so transpose.
        init_bw <-
            sample_size ^ {
                -2 / (total_num_vars + 4)
            } * sample_cov_hat
        init_bw_chol <- t(chol(init_bw))

        ## Set Cholesky decomposition and vectorized Cholesky decomposition in theta
        new_theta$bw_chol_decomp <- init_bw_chol

        lower_triangular_inds <- as.matrix(expand.grid(seq_len(ncol(
            new_theta$bw_chol_decomp
        )),
        seq_len(ncol(
            new_theta$bw_chol_decomp
        ))))
        lower_triangular_inds <-
            lower_triangular_inds[lower_triangular_inds[, 1] >= lower_triangular_inds[, 2], , drop = FALSE]

        new_theta$bw_chol_decomp_vec <-
            new_theta$bw_chol_decomp[lower_triangular_inds]

        ## Compute other bandwidth parameters used in calls to pdtmvn kernel
        bw_params <-
            compute_pdtmvn_kernel_bw_params_from_bw_chol_decomp(
                bw_chol_decomp_vec = new_theta$bw_chol_decomp_vec,
                bw_chol_decomp = new_theta$bw_chol_decomp,
                continuous_discrete_var_col_inds$continuous_vars,
                continuous_discrete_var_col_inds$discrete_vars
            )

        new_theta[names(bw_params)] <- bw_params

        new_theta$continuous_var_col_inds = continuous_discrete_var_col_inds$continuous_vars
        new_theta$discrete_var_col_inds = continuous_discrete_var_col_inds$discrete_vars
    } else {
        stop("Invalid parameterization for pdtmvn kernel function")
    }

    return(new_theta)
}
