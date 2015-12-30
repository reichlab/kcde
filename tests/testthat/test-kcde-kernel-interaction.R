library(kcde)
library(magrittr)
library(plyr)
library(dplyr)
library(mvtnorm)

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


test_that("initialize_theta works", {
        test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
        vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
            lag_value = c(1, 0, 2, 3, 2),
            stringsAsFactors = FALSE)
        vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
        
        leading_rows_to_drop <- 3
        trailing_rows_to_drop <- 1
        
        train_lagged_obs <- compute_lagged_obs_vecs(test_data,
            vars_and_lags,
            c(seq_len(leading_rows_to_drop),
                seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                    to=nrow(test_data))))
        
        expected <- list(
            initialize_params_pdtmvn_kernel(
                parameterization = "bw-diagonalized-est-eigenvalues",
                x = train_lagged_obs[, 1:2],
                continuous_vars = c("a_lag1", "b_lag0"),
                discrete_vars = NULL,
                discrete_var_range_fns = NULL,
                lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
                upper = c(a_lag1 = Inf, b_lag0 = Inf)),
            initialize_params_pdtmvn_kernel(
                parameterization = "bw-diagonalized-est-eigenvalues",
                x = train_lagged_obs[, 3:5],
                continuous_vars = c("b_lag2", "b_lag3"),
                discrete_vars = "c_lag2",
                discrete_var_range_fns = list(
                    c_lag2 = list(a = pdtmvn::floor_x_minus_1,
                        b = floor,
                        in_range = pdtmvn::equals_integer)),
                lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
                upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf))
        )
        
        
        kernel_components <- list(list(
                vars_and_lags = vars_and_lags[1:2, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[1:2],
                    discrete_vars = NULL,
                    discrete_var_range_fns = NULL,
                    lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
                    upper = c(a_lag1 = Inf, b_lag0 = Inf)
                ),
                vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
                vectorize_theta_est_args = NULL,
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ),
            list(
                vars_and_lags = vars_and_lags[3:5, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[3:4],
                    discrete_vars = vars_and_lags$combined_name[5],
                    discrete_var_range_fns = list(
                        c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
                    upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
                ),
                vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
                vectorize_theta_est_args = NULL,
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ))
        
        
        actual <- initialize_theta(prev_theta = NULL,
            updated_vars_and_lags = vars_and_lags,
            update_var_name = "b",
            update_lag_value = 3L,
            data = train_lagged_obs,
            kcde_control = list(kernel_components = kernel_components))
        
        expect_identical(actual, expected)
    })

test_that("extract_vectorized_theta_est_from_theta works", {
        vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
            lag_value = c(1, 0, 2, 3, 2),
            stringsAsFactors = FALSE)
        vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
        bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
        
        bw_eigen_1 <- eigen(diag(bws[1:2]))
        bw_eigen_2 <- eigen(diag(bws[3:5]))
        theta <- list(
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
                    bw_evals = bw_eigen_1$values,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = NULL,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL,
                    discrete_var_range_fns = NULL,
                    lower = rep(-Inf, 2),
                    upper = rep(Inf, 2),
                    log = TRUE
                )
            ),
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
                    bw_evals = bw_eigen_2$values,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = "c",
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3,
                    discrete_var_range_fns = list(
                        c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = rep(-Inf, 3),
                    upper = rep(Inf, 3),
                    log = TRUE
                )
            )
        )
        
        kernel_components <- list(list(
                vars_and_lags = vars_and_lags[1:2, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[1:2],
                    discrete_vars = NULL,
                    discrete_var_range_fns = NULL,
                    lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
                    upper = c(a_lag1 = Inf, b_lag0 = Inf)
                ),
                vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
                vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ),
            list(
                vars_and_lags = vars_and_lags[3:5, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[3:4],
                    discrete_vars = vars_and_lags$combined_name[5],
                    discrete_var_range_fns = list(
                        c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
                    upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
                ),
                vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
                vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ))
        
        expected <- c(log(bw_eigen_1$values), log(bw_eigen_2$values))
        
        actual <- extract_vectorized_theta_est_from_theta(theta = theta,
            vars_and_lags = vars_and_lags,
            kcde_control = list(kernel_components = kernel_components),
            add_fixed_params = FALSE)
        
        expect_identical(actual, expected)
    })

test_that("update_theta_from_vectorized_theta_est works", {
        vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
            lag_value = c(1, 0, 2, 3, 2),
            stringsAsFactors = FALSE)
        vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
        bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
        
        bw_eigen_1 <- eigen(diag(bws[1:2]))
        bw_eigen_2 <- eigen(diag(bws[3:5]))
        theta <- list(
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
                    bw_evals = bw_eigen_1$values,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = NULL,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL,
                    discrete_var_range_fns = NULL,
                    lower = rep(-Inf, 2),
                    upper = rep(Inf, 2),
                    log = TRUE
                )
            ),
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
                    bw_evals = bw_eigen_2$values,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = "c",
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3,
                    discrete_var_range_fns = list(
                        c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = rep(-Inf, 3),
                    upper = rep(Inf, 3),
                    log = TRUE
                )
            )
        )
        
        theta_est_vector <- (-11:-15)
        expected <- list(
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
                    bw_evals = exp(theta_est_vector[1:2]),
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = NULL,
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = NULL,
                    discrete_var_range_fns = NULL,
                    lower = rep(-Inf, 2),
                    upper = rep(Inf, 2),
                    log = TRUE
                )
            ),
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
                    bw_evals = exp(theta_est_vector[3:5]),
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = "c",
                    continuous_var_col_inds = 1:2,
                    discrete_var_col_inds = 3,
                    discrete_var_range_fns = list(
                        c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = rep(-Inf, 3),
                    upper = rep(Inf, 3),
                    log = TRUE
                )
            )
        )
        
        kernel_components <- list(list(
                vars_and_lags = vars_and_lags[1:2, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[1:2],
                    discrete_vars = NULL,
                    discrete_var_range_fns = NULL,
                    lower = c(a_lag1 = -Inf, b_lag0 = -Inf),
                    upper = c(a_lag1 = Inf, b_lag0 = Inf)
                ),
                vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
                vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ),
            list(
                vars_and_lags = vars_and_lags[3:5, ],
                kernel_fn = pdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues",
                    continuous_vars = vars_and_lags$combined_name[3:4],
                    discrete_vars = vars_and_lags$combined_name[5],
                    discrete_var_range_fns = list(
                        c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = c(b_lag2 = -Inf, b_lag3 = -Inf, c_lag2 = -Inf),
                    upper = c(b_lag2 = Inf, b_lag3 = Inf, c_lag2 = Inf)
                ),
                vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
                vectorize_kernel_params_args = list(parameterization = "bw-diagonalized-est-eigenvalues"),
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ))
        
        actual <- update_theta_from_vectorized_theta_est(theta_est_vector,
            theta,
            kcde_control = list(kernel_components = kernel_components))
        
        expect_identical(actual, expected)
    }
)

test_that("compute_kernel_values works -- one component, 1 continuous var, 1 discrete var, 1 var used", {
        test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
        vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
            lag_value = c(1, 0, 2, 3, 2),
            stringsAsFactors = FALSE)
        vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
        
        leading_rows_to_drop <- 3
        trailing_rows_to_drop <- 1
        
        train_lagged_obs <- compute_lagged_obs_vecs(test_data,
            vars_and_lags,
            c(seq_len(leading_rows_to_drop),
                seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                    to=nrow(test_data))))
        
        prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
            vars_and_lags,
            seq_len(9))
        
        bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
        expected <- pdtmvn::dpdtmvn(x = train_lagged_obs[, 2:5, drop = FALSE],
            mean = unlist(prediction_lagged_obs[, 2:5, drop = FALSE]),
            sigma = diag(bws[2:5]),
            continuous_vars = 1:3,
            discrete_vars = 4,
            log = TRUE)
        names(expected) <- NULL
        
        
        
        kernel_components <- list(list(
                vars_and_lags = vars_and_lags[2:5, ],
                kernel_fn = pdtmvn_kernel,
                rkernel_fn = rpdtmvn_kernel,
                theta_fixed = NULL,
                theta_est = list("bw"),
                initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
                initialize_kernel_params_args = list(
                    continuous_vars = vars_and_lags$combined_name[2:4],
                    discrete_vars = vars_and_lags$combined_name[5],
                    discrete_var_range_fns = list(
                        c = list(a = "pdtmvn::floor_x_minus_1", b = "floor", in_range = "pdtmvn::equals_integer"))
                ),
                vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
                vectorize_theta_est_args = NULL,
                update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
                update_theta_from_vectorized_theta_est_args = list(
                    parameterization = "bw-diagonalized-est-eigenvalues"
                )
            ))
        
        bw_eigen <- eigen(diag(bws[2:5]))
        theta <- list(
            c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen$vectors,
                    bw_evals = bw_eigen$values,
                    continuous_var_col_inds = 1:3,
                    discrete_var_col_inds = 4),
                list(
                    continuous_vars = vars_and_lags$combined_name,
                    discrete_vars = "c",
                    continuous_var_col_inds = 1:3,
                    discrete_var_col_inds = 4,
                    discrete_var_range_fns = list(
                        c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                    lower = rep(-Inf, 4),
                    upper = rep(Inf, 4),
                    log = TRUE
                )
            )
        )
        
        actual <- compute_kernel_values(train_obs = train_lagged_obs,
            prediction_obs = prediction_lagged_obs,
            kernel_components = kernel_components,
            theta = theta,
            log = TRUE)
        
        expect_identical(actual, expected)
    }
)


test_that("compute_kernel_values works -- two components", {
    test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
    vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
        lag_value = c(1, 0, 2, 3, 2),
        stringsAsFactors = FALSE)
    vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
    
    leading_rows_to_drop <- 3
    trailing_rows_to_drop <- 1
    
    train_lagged_obs <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        seq_len(9))
    
    bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
    expected <-
        pdtmvn::dpdtmvn(x = train_lagged_obs[, 1:2, drop = FALSE],
            mean = unlist(prediction_lagged_obs[, 1:2, drop = FALSE]),
            sigma = diag(bws[1:2]),
            continuous_vars = 1:2,
            discrete_vars = NULL,
            log = TRUE) +
        pdtmvn::dpdtmvn(x = train_lagged_obs[, 3:5, drop = FALSE],
            mean = unlist(prediction_lagged_obs[, 3:5, drop = FALSE]),
            sigma = diag(bws[3:5]),
            continuous_vars = 1:2,
            discrete_vars = 3,
            log = TRUE)
    names(expected) <- NULL
    
    
    
    kernel_components <- list(list(
            vars_and_lags = vars_and_lags[1:2, ],
            kernel_fn = pdtmvn_kernel,
            rkernel_fn = rpdtmvn_kernel,
            theta_fixed = NULL,
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
            initialize_kernel_params_args = list(
                continuous_vars = vars_and_lags$combined_name[1:2],
                discrete_vars = NULL,
                discrete_var_range_fns = NULL
            ),
            vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
            vectorize_theta_est_args = NULL,
            update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
            update_theta_from_vectorized_theta_est_args = list(
                parameterization = "bw-diagonalized-est-eigenvalues"
            )
        ),
        list(
            vars_and_lags = vars_and_lags[3:5, ],
            kernel_fn = pdtmvn_kernel,
            rkernel_fn = rpdtmvn_kernel,
            theta_fixed = NULL,
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
            initialize_kernel_params_args = list(
                continuous_vars = vars_and_lags$combined_name[3:4],
                discrete_vars = vars_and_lags$combined_name[5],
                discrete_var_range_fns = list(
                    c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer))
            ),
            vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
            vectorize_theta_est_args = NULL,
            update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
            update_theta_from_vectorized_theta_est_args = list(
                parameterization = "bw-diagonalized-est-eigenvalues"
            )
        ))
    
    bw_eigen_1 <- eigen(diag(bws[1:2]))
    bw_eigen_2 <- eigen(diag(bws[3:5]))
    theta <- list(
        c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
                bw_evals = bw_eigen_1$values,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = NULL),
            list(
                continuous_vars = vars_and_lags$combined_name,
                discrete_vars = NULL,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = NULL,
                discrete_var_range_fns = NULL,
                lower = rep(-Inf, 2),
                upper = rep(Inf, 2),
                log = TRUE
            )
        ),
        c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
                bw_evals = bw_eigen_2$values,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = 3),
            list(
                continuous_vars = vars_and_lags$combined_name,
                discrete_vars = "c",
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = 3,
                discrete_var_range_fns = list(
                    c = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer)),
                lower = rep(-Inf, 3),
                upper = rep(Inf, 3),
                log = TRUE
            )
        )
    )
    
    actual <- compute_kernel_values(train_obs = train_lagged_obs,
        prediction_obs = prediction_lagged_obs,
        kernel_components = kernel_components,
        theta = theta,
        log = TRUE)
    
    expect_identical(actual, expected)
})



test_that("simulate_values_from_product_kernel works -- two components", {
    ## Used as discretization function --
    ## consistently rounds x so that _.5 -> _+1
    round_up_.5 <- function(x) {
        return(sapply(x, function(x_i) {
                    if(x_i - floor(x_i) >= 0.5) {
                        return(ceiling(x_i))
                    } else {
                        return(floor(x_i))
                    }
                }))
    }
    
    test_data <- data.frame(a = rnorm(10), b = rnorm(10), c = 1:10)
    vars_and_lags <- data.frame(var_name = c("a", "b", "b", "b", "c"),
        lag_value = c(1, 0, 2, 3, 2),
        stringsAsFactors = FALSE)
    vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
    
    leading_rows_to_drop <- 3
    trailing_rows_to_drop <- 1
    
    train_lagged_obs <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        seq_len(9))
    
    bws <- c(1.3, 1.4, 1.5, 0.8, 0.34)
    n <- 10
    set.seed(1234)
    
    expected <- cbind(pdtmvn::rpdtmvn(n = n,
            x_fixed = prediction_lagged_obs[, 1, drop = FALSE],
            mean = train_lagged_obs[1, 1:2, drop = FALSE],
            sigma = diag(bws[1:2]),
            continuous_vars = 1:2,
            discrete_vars = NULL,
            validate_level = 1L),
        pdtmvn::rpdtmvn(n = n,
            x_fixed = prediction_lagged_obs[, 3, drop = FALSE],
            mean = train_lagged_obs[1, 3:5, drop = FALSE],
            sigma = diag(bws[3:5]),
            continuous_vars = 1:2,
            discrete_vars = 3,
            discrete_var_range_fns = list(
                c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5)),
            validate_level = 1L)
    )
    
    kernel_components <- list(list(
            vars_and_lags = vars_and_lags[1:2, ],
            kernel_fn = pdtmvn_kernel,
            rkernel_fn = rpdtmvn_kernel,
            theta_fixed = NULL,
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
            initialize_kernel_params_args = list(
                continuous_vars = vars_and_lags$combined_name[1:2],
                discrete_vars = NULL,
                discrete_var_range_fns = NULL
            ),
            vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
            vectorize_theta_est_args = NULL,
            update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
            update_theta_from_vectorized_theta_est_args = list(
                parameterization = "bw-diagonalized-est-eigenvalues"
            )
        ),
        list(
            vars_and_lags = vars_and_lags[3:5, ],
            kernel_fn = pdtmvn_kernel,
            rkernel_fn = rpdtmvn_kernel,
            theta_fixed = NULL,
            theta_est = list("bw"),
            initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
            initialize_kernel_params_args = list(
                continuous_vars = vars_and_lags$combined_name[3:4],
                discrete_vars = vars_and_lags$combined_name[5],
                discrete_var_range_fns = list(
                    c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5))
            ),
            vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
            vectorize_theta_est_args = NULL,
            update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
            update_theta_from_vectorized_theta_est_args = list(
                parameterization = "bw-diagonalized-est-eigenvalues"
            )
        ))
    
    bw_eigen_1 <- eigen(diag(bws[1:2]))
    bw_eigen_2 <- eigen(diag(bws[3:5]))
    theta <- list(
        c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_1$vectors,
                bw_evals = bw_eigen_1$values,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = NULL),
            list(
                continuous_vars = vars_and_lags$combined_name[1:2],
                discrete_vars = NULL,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = NULL,
                discrete_var_range_fns = NULL,
                lower = rep(-Inf, 2),
                upper = rep(Inf, 2),
                log = TRUE
            )
        ),
        c(compute_pdtmvn_kernel_bw_params_from_bw_eigen(bw_evecs = bw_eigen_2$vectors,
                bw_evals = bw_eigen_2$values,
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = 3),
            list(
                continuous_vars = vars_and_lags$combined_name[3:4],
                discrete_vars = "c_lag2",
                continuous_var_col_inds = 1:2,
                discrete_var_col_inds = 3,
                discrete_var_range_fns = list(
                    c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5)),
                lower = rep(-Inf, 3),
                upper = rep(Inf, 3),
                log = TRUE
            )
        )
    )
    
    set.seed(1234)
    actual <- simulate_values_from_product_kernel(n = n,
        conditioning_obs = prediction_lagged_obs[, c(1, 3), drop = FALSE],
        center = train_lagged_obs[1, , drop = FALSE],
        kernel_components = kernel_components,
        theta = theta)
    
    expect_identical(actual, expected)
})
