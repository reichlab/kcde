## all functions in kcde-kernel-interaction.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
##
## X initialize_theta
## X extract_vectorized_theta_est_from_theta
## X update_theta_from_vectorized_theta_est
## X compute_kernel_values
##   simulate_values_from_product_kernel


context("kcde-kernel-interaction")


## test_that("initialize_theta works", {
##         test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
##         vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##             offset_value = c(1, 0, 2, 3, 2),
##             offset_type = rep("lag", 5),
##             stringsAsFactors = FALSE)
##         vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_lag", vars_and_offsets$offset_value)
        
##         leading_rows_to_drop <- 3
##         trailing_rows_to_drop <- 1
        
##         train_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##             vars_and_offsets = vars_and_offsets,
##             leading_rows_to_drop = leading_rows_to_drop,
##             trailing_rows_to_drop = trailing_rows_to_drop)
        
##         expected_2 <- list(
##             initialize_params_pdtmvn_kernel(
##                 prev_theta = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = c("a_lag1", "b_lag0"),
##                     discrete_vars = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
##                     upper = c(a_lag1 = Inf, b_lag0 = Inf)),
##                 x = train_lagged_obs[, 1:2]),
##             initialize_params_pdtmvn_kernel(
##                 prev_theta = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = c("b_lag2", "b_lag3"),
##                     discrete_vars = "c_lag2",
##                     discrete_var_range_fns = list(
##                         c_lag2 = list(a = pdtmvn::floor_x_minus_1,
##                             b = floor,
##                             in_range = pdtmvn::equals_integer)),
##                     lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
##                     upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)),
##                 x = train_lagged_obs[, 3:5])
##         )
        
##         expected_1 <- expected_2
##         expected_1[1] <- list(NULL)
        
        
##         kernel_components <- list(list(
##                 vars_and_offsets = vars_and_offsets[1:2, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[1:2],
##                     discrete_vars = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
##                     upper = c(a_lag1 = Inf, b_lag0 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_theta_est_args = NULL,
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ),
##             list(
##                 vars_and_offsets = vars_and_offsets[3:5, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[3:4],
##                     discrete_vars = vars_and_offsets$combined_name[5],
##                     discrete_var_range_fns = list(
##                         c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
##                     upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_theta_est_args = NULL,
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ))
        
        
##         actual_1 <- initialize_theta(prev_theta = NULL,
##             updated_vars_and_offsets = vars_and_offsets,
##             update_var_name = "b",
##             update_offset_value = 3L,
##             update_offset_type = "lag",
##             data = train_lagged_obs,
##             kcde_control = list(kernel_components = kernel_components))
        
##         expect_identical(actual_1, expected_1)
        
        
##         actual_2 <- initialize_theta(prev_theta = actual_1,
##             updated_vars_and_offsets = vars_and_offsets,
##             update_var_name = "b",
##             update_offset_value = 0L,
##             update_offset_type = "lag",
##             data = train_lagged_obs,
##             kcde_control = list(kernel_components = kernel_components))
        
##         expect_identical(actual_2, expected_2)
##     })

## test_that("extract_vectorized_theta_est_from_theta works", {
##         vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##             offset_value = c(1, 0, 2, 3, 2),
##             offset_type = rep("lag", 5),
##             stringsAsFactors = FALSE)
##         vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_", vars_and_offsets$offset_type, vars_and_offsets$offset_value)
##         bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
        
##         bw_eigen_1 <- eigen(diag(bws[1:2]))
##         bw_eigen_2 <- eigen(diag(bws[3:5]))
##         theta <- list(
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
##                     bw_evals = bw_eigen_1$values,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = NULL,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = rep(-Inf, 2),
##                     upper = rep(Inf, 2),
##                     log = TRUE
##                 )
##             ),
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
##                     bw_evals = bw_eigen_2$values,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = "c",
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3,
##                     discrete_var_range_fns = list(
##                         c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = rep(-Inf, 3),
##                     upper = rep(Inf, 3),
##                     log = TRUE
##                 )
##             )
##         )
        
##         kernel_components <- list(list(
##                 vars_and_offsets = vars_and_offsets[1:2, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[1:2],
##                     discrete_vars = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
##                     upper = c(a_lag1 = Inf, b_lag0 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ),
##             list(
##                 vars_and_offsets = vars_and_offsets[3:5, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[3:4],
##                     discrete_vars = vars_and_offsets$combined_name[5],
##                     discrete_var_range_fns = list(
##                         c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
##                     upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ))
        
##         expected <- c(log(bw_eigen_1$values), log(bw_eigen_2$values))
        
##         actual <- extract_vectorized_theta_est_from_theta(theta = theta,
##             vars_and_offsets = vars_and_offsets,
##             kcde_control = list(kernel_components = kernel_components))
        
##         expect_identical(actual, expected)
##     })

## test_that("update_theta_from_vectorized_theta_est works", {
##         vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##             offset_value = c(1, 0, 2, 3, 2),
##             offset_type = rep("lag", 5),
##             stringsAsFactors = FALSE)
##         vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_", vars_and_offsets$offset_type, vars_and_offsets$offset_value)
##         bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
        
##         bw_eigen_1 <- eigen(diag(bws[1:2]))
##         bw_eigen_2 <- eigen(diag(bws[3:5]))
##         theta <- list(
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
##                     bw_evals = bw_eigen_1$values,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = NULL,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = rep(-Inf, 2),
##                     upper = rep(Inf, 2),
##                     log = TRUE
##                 )
##             ),
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
##                     bw_evals = bw_eigen_2$values,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = "c",
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3,
##                     discrete_var_range_fns = list(
##                         c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = rep(-Inf, 3),
##                     upper = rep(Inf, 3),
##                     log = TRUE
##                 )
##             )
##         )
        
##         theta_est_vector <- (-11:-15)
##         expected <- list(
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
##                     bw_evals = exp(theta_est_vector[1:2]),
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = NULL,
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = rep(-Inf, 2),
##                     upper = rep(Inf, 2),
##                     log = TRUE
##                 )
##             ),
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
##                     bw_evals = exp(theta_est_vector[3:5]),
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = "c",
##                     continuous_var_col_inds = 1:2,
##                     discrete_var_col_inds = 3,
##                     discrete_var_range_fns = list(
##                         c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = rep(-Inf, 3),
##                     upper = rep(Inf, 3),
##                     log = TRUE
##                 )
##             )
##         )
        
##         kernel_components <- list(list(
##                 vars_and_offsets = vars_and_offsets[1:2, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[1:2],
##                     discrete_vars = NULL,
##                     discrete_var_range_fns = NULL,
##                     lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
##                     upper = c(a_lag1 = Inf, b_lag0 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ),
##             list(
##                 vars_and_offsets = vars_and_offsets[3:5, ],
##                 kernel_fn = pdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[3:4],
##                     discrete_vars = vars_and_offsets$combined_name[5],
##                     discrete_var_range_fns = list(
##                         c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
##                     upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ))
        
##         actual <- update_theta_from_vectorized_theta_est(theta_est_vector,
##             theta,
##             kcde_control = list(kernel_components = kernel_components))
        
##         expect_identical(actual, expected)
##     }
## )

## test_that("compute_kernel_values works -- one component, 1 continuous var, 1 discrete var, 1 var used", {
##         test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
##         vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##             offset_value = c(1, 0, 2, 3, 2),
##             offset_type = rep("lag", 5),
##             stringsAsFactors = FALSE)
##         vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_", vars_and_offsets$offset_type, vars_and_offsets$offset_value)
        
##         leading_rows_to_drop <- 3
##         trailing_rows_to_drop <- 1
        
##         train_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##             vars_and_offsets = vars_and_offsets,
##             leading_rows_to_drop = leading_rows_to_drop,
##             trailing_rows_to_drop = trailing_rows_to_drop)
        
##         prediction_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##             vars_and_offsets = vars_and_offsets,
##             leading_rows_to_drop = 9,
##             trailing_rows_to_drop = 0)
        
##         bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
##         expected <- pdtmvn::dpdtmvn(x = train_lagged_obs[, 2:5, drop = FALSE],
##             mean = unlist(prediction_lagged_obs[, 2:5, drop = FALSE]),
##             sigma = diag(bws[2:5]),
##             continuous_vars = 1:3,
##             discrete_vars = 4,
##             log = TRUE)
##         names(expected) <- NULL
        
        
        
##         kernel_components <- list(list(
##                 vars_and_offsets = vars_and_offsets[2:5, ],
##                 kernel_fn = pdtmvn_kernel,
##                 rkernel_fn = rpdtmvn_kernel,
##                 theta_fixed = list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name[2:4],
##                     discrete_vars = vars_and_offsets$combined_name[5],
##                     discrete_var_range_fns = list(
##                         c = list(a = "pdtmvn::floor_x_minus_1", b = "floor", in_range = "pdtmvn::equals_integer"))
##                 ),
##                 theta_est = list("bw"),
##                 initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##                 initialize_kernel_params_args = NULL,
##                 vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##                 vectorize_theta_est_args = NULL,
##                 update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##                 update_theta_from_vectorized_theta_est_args = NULL
##             ))
        
##         bw_eigen <- eigen(diag(bws[2:5]))
##         theta <- list(
##             c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen$vectors,
##                     bw_evals = bw_eigen$values,
##                     continuous_var_col_inds = 1:3,
##                     discrete_var_col_inds = 4),
##                 list(
##                     parameterization = "bw-diagonalized-est-eigenvalues",
##                     continuous_vars = vars_and_offsets$combined_name,
##                     discrete_vars = "c",
##                     continuous_var_col_inds = 1:3,
##                     discrete_var_col_inds = 4,
##                     discrete_var_range_fns = list(
##                         c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                     lower = rep(-Inf, 4),
##                     upper = rep(Inf, 4),
##                     log = TRUE
##                 )
##             )
##         )
        
##         actual <- compute_kernel_values(train_obs = train_lagged_obs,
##             prediction_obs = prediction_lagged_obs,
##             kernel_components = kernel_components,
##             theta = theta,
##             log = TRUE)
        
##         expect_identical(actual, expected)
##     }
## )


## test_that("compute_kernel_values works -- two components", {
##     test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
##     vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##         offset_value = c(1, 0, 2, 3, 2),
##         offset_type = rep("lag", 5),
##         stringsAsFactors = FALSE)
##     vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_", vars_and_offsets$offset_type, vars_and_offsets$offset_value)
    
##     leading_rows_to_drop <- 3
##     trailing_rows_to_drop <- 1
    
##     train_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##         vars_and_offsets = vars_and_offsets,
##         leading_rows_to_drop = leading_rows_to_drop,
##         trailing_rows_to_drop = trailing_rows_to_drop)
    
##     prediction_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##         vars_and_offsets = vars_and_offsets,
##         leading_rows_to_drop = 9,
##         trailing_rows_to_drop = 0)
    
##     bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
##     expected <-
##         pdtmvn::dpdtmvn(x = train_lagged_obs[, 1:2, drop = FALSE],
##             mean = unlist(prediction_lagged_obs[, 1:2, drop = FALSE]),
##             sigma = diag(bws[1:2]),
##             continuous_vars = 1:2,
##             discrete_vars = NULL,
##             log = TRUE) +
##         pdtmvn::dpdtmvn(x = train_lagged_obs[, 3:5, drop = FALSE],
##             mean = unlist(prediction_lagged_obs[, 3:5, drop = FALSE]),
##             sigma = diag(bws[3:5]),
##             continuous_vars = 1:2,
##             discrete_vars = 3,
##             log = TRUE)
##     names(expected) <- NULL
    
    
    
##     kernel_components <- list(list(
##             vars_and_offsets = vars_and_offsets[1:2, ],
##             kernel_fn = pdtmvn_kernel,
##             rkernel_fn = rpdtmvn_kernel,
##             theta_fixed = list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name[1:2],
##                 discrete_vars = NULL,
##                 discrete_var_range_fns = NULL
##             ),
##             theta_est = list("bw"),
##             initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##             initialize_kernel_params_args = NULL,
##             vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##             vectorize_theta_est_args = NULL,
##             update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##             update_theta_from_vectorized_theta_est_args = NULL
##         ),
##         list(
##             vars_and_offsets = vars_and_offsets[3:5, ],
##             kernel_fn = pdtmvn_kernel,
##             rkernel_fn = rpdtmvn_kernel,
##             theta_fixed = list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name[3:4],
##                 discrete_vars = vars_and_offsets$combined_name[5],
##                 discrete_var_range_fns = list(
##                     c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer))
##             ),
##             theta_est = list("bw"),
##             initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##             initialize_kernel_params_args = NULL,
##             vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##             vectorize_theta_est_args = NULL,
##             update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##             update_theta_from_vectorized_theta_est_args = NULL
##         ))
    
##     bw_eigen_1 <- eigen(diag(bws[1:2]))
##     bw_eigen_2 <- eigen(diag(bws[3:5]))
##     theta <- list(
##         c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
##                 bw_evals = bw_eigen_1$values,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = NULL),
##             list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name,
##                 discrete_vars = NULL,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = NULL,
##                 discrete_var_range_fns = NULL,
##                 lower = rep(-Inf, 2),
##                 upper = rep(Inf, 2),
##                 log = TRUE
##             )
##         ),
##         c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
##                 bw_evals = bw_eigen_2$values,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = 3),
##             list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name,
##                 discrete_vars = "c",
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = 3,
##                 discrete_var_range_fns = list(
##                     c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
##                 lower = rep(-Inf, 3),
##                 upper = rep(Inf, 3),
##                 log = TRUE
##             )
##         )
##     )
    
##     actual <- compute_kernel_values(train_obs = train_lagged_obs,
##         prediction_obs = prediction_lagged_obs,
##         kernel_components = kernel_components,
##         theta = theta,
##         log = TRUE)
    
##     expect_identical(actual, expected)
## })



## test_that("simulate_values_from_product_kernel works -- two components", {
##     ## Used as discretization function --
##     ## consistently rounds x so that _.5 -> _+1
##     round_up_.5 <- function(x) {
##         return(sapply(x, function(x_i) {
##                     if(x_i - floor(x_i) >= 0.5) {
##                         return(ceiling(x_i))
##                     } else {
##                         return(floor(x_i))
##                     }
##                 }))
##     }
    
##     test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
##     vars_and_offsets <- data.frame(var_name = c("a", "b", "b", "b", "c"),
##         offset_value = c(1, 0, 2, 3, 2),
##         offset_type = rep("lag", 5),
##         stringsAsFactors = FALSE)
##     vars_and_offsets$combined_name <- paste0(vars_and_offsets$var_name, "_", vars_and_offsets$offset_type, vars_and_offsets$offset_value)
    
##     leading_rows_to_drop <- 3
##     trailing_rows_to_drop <- 1
    
##     train_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##         vars_and_offsets = vars_and_offsets,
##         leading_rows_to_drop = leading_rows_to_drop,
##         trailing_rows_to_drop = trailing_rows_to_drop)
    
##     prediction_lagged_obs <- compute_offset_obs_vecs(data = test_data,
##         vars_and_offsets = vars_and_offsets,
##         leading_rows_to_drop = 9,
##         trailing_rows_to_drop = 0)
    
##     bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
##     n <- 10
##     set.seed(1234)
    
##     expected <- cbind(pdtmvn::rpdtmvn(n = n,
##             x_fixed = prediction_lagged_obs[, 1, drop = FALSE],
##             mean = train_lagged_obs[1, 1:2, drop = FALSE],
##             sigma = diag(bws[1:2]),
##             continuous_vars = 1:2,
##             discrete_vars = NULL,
##             validate_level = 1L),
##         pdtmvn::rpdtmvn(n = n,
##             x_fixed = prediction_lagged_obs[, 3, drop = FALSE],
##             mean = train_lagged_obs[1, 3:5, drop = FALSE],
##             sigma = diag(bws[3:5]),
##             continuous_vars = 1:2,
##             discrete_vars = 3,
##             discrete_var_range_fns = list(
##                 c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5)),
##             validate_level = 1L)
##     )
    
##     kernel_components <- list(list(
##             vars_and_offsets = vars_and_offsets[1:2, ],
##             kernel_fn = pdtmvn_kernel,
##             rkernel_fn = rpdtmvn_kernel,
##             theta_fixed = list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name[1:2],
##                 discrete_vars = NULL,
##                 discrete_var_range_fns = NULL
##             ),
##             theta_est = list("bw"),
##             initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##             initialize_kernel_params_args = NULL,
##             vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##             vectorize_theta_est_args = NULL,
##             update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##             update_theta_from_vectorized_theta_est_args = NULL
##         ),
##         list(
##             vars_and_offsets = vars_and_offsets[3:5, ],
##             kernel_fn = pdtmvn_kernel,
##             rkernel_fn = rpdtmvn_kernel,
##             theta_fixed = list(
##                 parameterization = "bw-diagonalized-est-eigenvalues",
##                 continuous_vars = vars_and_offsets$combined_name[3:4],
##                 discrete_vars = vars_and_offsets$combined_name[5],
##                 discrete_var_range_fns = list(
##                     c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5))
##             ),
##             theta_est = list("bw"),
##             initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
##             initialize_kernel_params_args = NULL,
##             vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
##             vectorize_theta_est_args = NULL,
##             update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
##             update_theta_from_vectorized_theta_est_args = NULL
##         ))
    
##     bw_eigen_1 <- eigen(diag(bws[1:2]))
##     bw_eigen_2 <- eigen(diag(bws[3:5]))
##     theta <- list(
##         c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
##                 bw_evals = bw_eigen_1$values,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = NULL),
##             list(
##                 continuous_vars = vars_and_offsets$combined_name[1:2],
##                 discrete_vars = NULL,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = NULL,
##                 discrete_var_range_fns = NULL,
##                 lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
##                 upper = c(a_lag1 = Inf, b_lag0 = Inf),
##                 log = TRUE
##             )
##         ),
##         c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
##                 bw_evals = bw_eigen_2$values,
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = 3),
##             list(
##                 continuous_vars = vars_and_offsets$combined_name[3:4],
##                 discrete_vars = "c_lag2",
##                 continuous_var_col_inds = 1:2,
##                 discrete_var_col_inds = 3,
##                 discrete_var_range_fns = list(
##                     c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5)),
##                 lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
##                 upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf),
##                 log = TRUE
##             )
##         )
##     )
    
##     set.seed(1234)
##     actual <- simulate_values_from_product_kernel(n = n,
##         conditioning_obs = prediction_lagged_obs[, c(1, 3), drop = FALSE],
##         center = train_lagged_obs[1, , drop = FALSE],
##         kernel_components = kernel_components,
##         theta = theta)
    
##     expect_identical(actual, expected)
## })
