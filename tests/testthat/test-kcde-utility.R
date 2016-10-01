## all functions in kcde-utility.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
##
## X update_vars_and_lags
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

test_that("compute_offset_obs_vecs works -- all variables used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    vars_and_offsets <- data.frame(var_name = c("a", "a", "b", "b"),
                                   offset_value = c(1, 5, 0, 2),
                                   offset_type = c("horizon", "lag", "lag", "lag"),
                                   stringsAsFactors = FALSE)
    vars_and_offsets$combined_name <-
        paste0(vars_and_offsets$var_name,
               "_",
               vars_and_offsets$offset_type,
               vars_and_offsets$offset_value)
    ## lags <- list(b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 1

    expected <- test_data
    expected <- expected %>% mutate(a_horizon1 = lead(a, 1),
                                    a_lag5 = lag(a, 5),
                                    b_lag0 = b,
                                    b_lag2 = lag(b, 2))
    expected <- expected[seq(from = 7, to = 9), 3:6]

    actual <- compute_offset_obs_vecs(data = test_data,
                                      filter_control = NULL,
                                      phi = NULL,
                                      vars_and_offsets = vars_and_offsets,
                                      time_name = NULL,
                                      leading_rows_to_drop = leading_rows_to_drop,
                                      trailing_rows_to_drop = trailing_rows_to_drop,
                                      additional_rows_to_drop = NULL,
                                      na.action = "na.omit")

    expect_identical(actual, expected)
})

test_that("compute_offset_obs_vecs works -- one variable not used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
	vars_and_offsets <- data.frame(var_name = c("b", "b"),
                                   offset_value = c(0, 2),
                                   offset_type = c("lag", "lag"),
                                   stringsAsFactors = FALSE)
    vars_and_offsets$combined_name <- 
        paste0(vars_and_offsets$var_name,
               "_",
               vars_and_offsets$offset_type,
               vars_and_offsets$offset_value)
    ## lags <- list(b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 0

    expected <- test_data
    expected <- expected %>% mutate(b_lag0 = b,
                                    b_lag2 = lag(b, 2))
    expected <- expected[seq(from = 7, to = 10), 3:4]

    actual <- compute_offset_obs_vecs(data = test_data,
                                      filter_control = NULL,
                                      phi = NULL,
                                      vars_and_offsets = vars_and_offsets,
                                      time_name = NULL,
                                      leading_rows_to_drop = leading_rows_to_drop,
                                      trailing_rows_to_drop = 0,
                                      additional_rows_to_drop = NULL,
                                      na.action = "na.omit")

    expect_identical(actual, expected)
})

test_that("update_vars_and_offsets works -- removal", {
    orig_vars_and_offsets <-  data.frame(var_name = c("a", "b", "b", "b", "c"),
                                         offset_value = as.integer(c(1, 0, 2, 3, 2)),
                                         offset_type = rep("lag", 5),
                                         stringsAsFactors = FALSE)
    orig_vars_and_offsets$combined_name <- 
        paste0(orig_vars_and_offsets$var_name,
               "_",
               orig_vars_and_offsets$offset_type,
               orig_vars_and_offsets$offset_value)

    update_var_name <- "b"
    update_offset_value <- 3
    update_offset_type <- "lag"

    expected <- orig_vars_and_offsets[-4, , drop = FALSE]
    rownames(expected) <- NULL

    actual <- update_vars_and_offsets(prev_vars_and_offsets = orig_vars_and_offsets,
                                      update_var_name = update_var_name,
                                      update_offset_value = update_offset_value,
                                      update_offset_type = update_offset_type)

    expect_identical(actual, expected)
})

test_that("update_vars_and_offsets works -- addition", {
    orig_vars_and_offsets <-  data.frame(var_name = c("a", "b", "b", "b", "c"),
                                         offset_value = as.integer(c(1, 0, 2, 3, 2)),
                                         offset_type = rep("lag", 5),
                                         stringsAsFactors = FALSE)
    orig_vars_and_offsets$combined_name <- 
        paste0(orig_vars_and_offsets$var_name,
               "_",
               orig_vars_and_offsets$offset_type,
               orig_vars_and_offsets$offset_value)

    update_var_name <- "a"
    update_offset_value <- 4
    update_offset_type <- "horizon"

    expected <- orig_vars_and_offsets
    expected[6, ] <- c("a", 4L, "horizon", "a_horizon4")
    expected$offset_value <- as.integer(expected$offset_value)

    row_reordering <- order(expected$combined_name)
    expected <- expected[row_reordering, , drop = FALSE]
    rownames(expected) <- NULL

    actual <- update_vars_and_offsets(prev_vars_and_offsets = orig_vars_and_offsets,
                                      update_var_name = update_var_name,
                                      update_offset_value = update_offset_value,
                                      update_offset_type = update_offset_type)

    expect_identical(actual, expected)
})
