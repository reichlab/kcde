library(kcde)
library(magrittr)
library(plyr)
library(dplyr)
library(mvtnorm)

## all functions in kcde-utility.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
##
## X update_vars_and_lags
## X assemble_training_examples
##   assemble_prediction_examples
## X compute_normalized_log_weights
## X compute_lagged_obs_vecs
##   mase_from_kernel_weights_and_centers
##   mae_from_kernel_weights_and_centers
##   get_inds_smallest_k
##   logspace_sub
##   logspace_add
##   logspace_sum
##   logspace_sum_matrix_rows
##   logspace_sub_matrix_rows
##   mase
##   mae



context("kcde-utility")

test_that("compute_normalized_log_weights works", {
    init_val <- c(0, -100, -50, -78, -99, -0.2)
    
    result <- compute_normalized_log_weights(init_val)
    
    expect_equal(sum(exp(result)), 1)
    expect_true(all.equal(init_val - result,
        rep(init_val[1] - result[1], length(init_val))))
})

test_that("assemble_training_examples works", {
    test_data <- data.frame(a = 1:20, b = rnorm(20))
	vars_and_lags <- data.frame(var_name = c("a", "a", "b", "b"),
		lag_value = c(1, 5, 0, 2),
		stringsAsFactors = FALSE)
	vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
	leading_rows_to_drop <- 6

    actual <- assemble_training_examples(test_data,
        vars_and_lags,
        y_names="b",
        leading_rows_to_drop=leading_rows_to_drop,
        additional_rows_to_drop=c(14, 15, 18),
        prediction_horizon = 2,
        drop_trailing_rows=TRUE)
    
    expected_lagged <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        c(seq_len(leading_rows_to_drop), c(14, 15, 18, 19, 20)))
    expected_lead <- test_data[1:20 + 2, "b", drop=FALSE]
    expected_lead <- expected_lead[
        -c(seq_len(leading_rows_to_drop), c(14, 15, 18, 19, 20)), , drop=FALSE]
	colnames(expected_lead) <- "b_horizon2"
	rownames(expected_lead) <- as.integer(rownames(expected_lead)) - 2L
    expected <- list(lagged_obs=expected_lagged,
        lead_obs=expected_lead
    )
    
    expect_identical(actual, expected)
})

test_that("compute_lagged_obs_vecs works -- all variables used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
	vars_and_lags <- data.frame(var_name = c("a", "a", "b", "b"),
		lag_value = c(1, 5, 0, 2),
		stringsAsFactors = FALSE)
	vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 1
    
    expected <- test_data
    expected <- expected %>% mutate(a_lag1=lag(a, 1),
        a_lag5=lag(a, 5),
        b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=9), 3:6]
    
    actual <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    expect_identical(actual, expected)
})

test_that("compute_lagged_obs_vecs works -- one variable not used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
	vars_and_lags <- data.frame(var_name = c("b", "b"),
		lag_value = c(0, 2),
		stringsAsFactors = FALSE)
	vars_and_lags$combined_name <- paste0(vars_and_lags$var_name, "_lag", vars_and_lags$lag_value)
#	lags <- list(b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 0
    
    expected <- test_data
    expected <- expected %>% mutate(b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=10), 3:4]
    
    actual <- compute_lagged_obs_vecs(test_data,
        vars_and_lags,
        seq_len(leading_rows_to_drop))
    
    expect_identical(actual, expected)
})



test_that("update_vars_and_lags works", {
    orig_vars_and_lags <-  data.frame(var_name = c("a", "b", "b", "b", "c"),
        lag_value = as.integer(c(1, 0, 2, 3, 2)),
        stringsAsFactors = FALSE)
    orig_vars_and_lags$combined_name <- paste0(orig_vars_and_lags$var_name, "_lag", orig_vars_and_lags$lag_value)
    
    ## test removal
    update_var_name <- "b"
    update_lag_value <- 3
    
    expected <- orig_vars_and_lags[-4, , drop = FALSE]
    
    actual <- update_vars_and_lags(orig_vars_and_lags,
        update_var_name,
        update_lag_value)
    
    expect_identical(actual, expected)
        
    ## test addition
    update_var_name <- "a"
    update_lag_value <- 4
    
    expected <- orig_vars_and_lags
    expected[6, ] <- c("a", 4L, "a_lag4")
    expected$lag_value <- as.integer(expected$lag_value)
    
    actual <- update_vars_and_lags(orig_vars_and_lags,
        update_var_name,
        update_lag_value)
    
    expect_identical(actual, expected)
})
