library(cdcfluview)
library(dplyr)
library(lubridate)
library(ggplot2)
library(grid)
library(kcde)
library(proftools)
library(doMC)


usflu<-get_flu_data("national", "ilinet", years=1997:2015)
ili_national <- transmute(usflu,
    region.type = REGION.TYPE,
    region = REGION,
    year = YEAR,
    week = WEEK,
    total_cases = as.numeric(X..WEIGHTED.ILI))
ili_national$time <- ymd(paste(ili_national$year, "01", "01", sep = "-"))
week(ili_national$time) <- ili_national$week
ili_national$time_index <- seq_len(nrow(ili_national))

str(ili_national)

## separate kernel components for lagged total cases and leading total cases
kernel_components <- list(
    list(
        vars_and_offsets = data.frame(var_name = "time_index",
            offset_value = 0L,
            offset_type = "lag",
            combined_name = "time_index_lag0",
            stringsAsFactors = FALSE),
        kernel_fn = periodic_kernel,
        theta_fixed = list(period=pi / 52.2),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_periodic_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_periodic_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_periodic_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ),
    list(
        vars_and_offsets = data.frame(var_name = c("total_cases", "total_cases"),
            offset_value = c(0L, 1L),
            offset_type = c("lag", "lag"),
            combined_name = c("total_cases_lag0", "total_cases_lag1"),
            stringsAsFactors = FALSE),
        kernel_fn = pdtmvn_kernel,
        rkernel_fn = rpdtmvn_kernel,
        theta_fixed = list(
            parameterization = "bw-diagonalized-est-eigenvalues",
            continuous_vars = c("total_cases_lag0", "total_cases_lag1"),
            discrete_vars = NULL,
            discrete_var_range_fns = NULL,
            lower = c(total_cases_lag0 = -Inf, total_cases_lag1 = -Inf),
            upper = c(total_cases_lag0 = Inf, total_cases_lag1 = Inf)
        ),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ),
    list(
        vars_and_offsets = data.frame(var_name = "total_cases",
            offset_value = 1L,
            offset_type = "horizon",
            combined_name = "total_cases_horizon1",
            stringsAsFactors = FALSE),
        kernel_fn = pdtmvn_kernel,
        rkernel_fn = rpdtmvn_kernel,
        theta_fixed = list(
            parameterization = "bw-diagonalized-est-eigenvalues",
            continuous_vars = "total_cases_horizon1",
            discrete_vars = NULL,
            discrete_var_range_fns = NULL,
            lower = c(total_cases_horizon1 = -Inf),
            upper = c(total_cases_horizon1 = Inf)
        ),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ))


kcde_control <- create_kcde_control(X_names = "time_index",
    y_names = "total_cases",
    time_name = "time",
    prediction_horizons = 1L,
    kernel_components = kernel_components,
    crossval_buffer = ymd("2010-01-01") - ymd("2009-01-01"),
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_type = "distribution",
    loss_args = NULL)

#tmp_file <- tempfile()
#Rprof(tmp_file, gc.profiling = TRUE, line.profiling = TRUE)

## set up parallelization
registerDoMC(cores=3)


## estimate parameters using only data up through 2014
flu_kcde_fit_orig_scale <- kcde(data = ili_national[ili_national$year <= 2014, ],
    kcde_control = kcde_control)

#Rprof(NULL)
#pd <- readProfileData(tmp_file)
#options(width = 300)
#hotPaths(pd, maxdepth = 27)


## sample from predictive distribution for first week in 2015
predictive_sample <- kcde_predict(kcde_fit = flu_kcde_fit_orig_scale,
        prediction_data = ili_national[ili_national$year == 2014 & ili_national$week == 53, , drop = FALSE],
        leading_rows_to_drop = 0,
        trailing_rows_to_drop = 1L,
        additional_training_rows_to_drop = NULL,
        prediction_type = "sample",
        n = 10000L)

ggplot() +
	geom_density(aes(x = predictive_sample)) +
	geom_vline(aes(xintercept = ili_national$total_cases[ili_national$year == 2015 & ili_national$week == 1]),
		colour = "red") +
	xlab("Total Cases") +
	ylab("Predictive Density") +
	ggtitle("Realized total cases vs. one week ahead predictive density\nWeek 1 of 2015") +
	theme_bw()



## obtain kernel weights and centers
debug(kcde_kernel_centers_and_weights_predict_given_lagged_obs)
predictive_kernel_weights_and_centers <- kcde_predict(kcde_fit = flu_kcde_fit_orig_scale,
    prediction_data = ili_national[ili_national$year == 2014 & ili_national$week == 53, , drop = FALSE],
    leading_rows_to_drop = 0,
    trailing_rows_to_drop = 1L,
    additional_training_rows_to_drop = NULL,
    prediction_type = "centers-and-weights")




ggplot() +
    geom_point(aes(x = predictive_kernel_weights_and_centers$centers[, 1],
        y = predictive_kernel_weights_and_centers$weights)) +
    theme_bw()


matching_predictive_value <- compute_offset_obs_vecs(data = flu_kcde_fit_orig_scale$train_data,
    vars_and_offsets = flu_kcde_fit_orig_scale$vars_and_offsets,
    time_name = flu_kcde_fit_orig_scale$kcde_control$time_name,
    leading_rows_to_drop = 0,
    trailing_rows_to_drop = 1L,
    additional_rows_to_drop = NULL,
    na.action = flu_kcde_fit_orig_scale$kcde_control$na.action)

ggplot() +
    geom_point(aes(y = predictive_kernel_weights_and_centers$centers[, 1],
            x = matching_predictive_value$total_cases_lag0,
            alpha = predictive_kernel_weights_and_centers$weights)) +
    theme_bw()







## One kernel component for lagged total cases and leading total cases
kernel_components <- list(
    list(
        vars_and_offsets = data.frame(var_name = "time_index",
            offset_value = 0L,
            offset_type = "lag",
            combined_name = "time_index_lag0",
            stringsAsFactors = FALSE),
        kernel_fn = periodic_kernel,
        theta_fixed = list(period=pi / 52.2),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_periodic_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_periodic_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_periodic_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ),
    list(
        vars_and_offsets = data.frame(var_name = c("total_cases", "total_cases", "total_cases"),
            offset_value = c(0L, 1L, 1L),
            offset_type = c("lag", "lag", "horizon"),
            combined_name = c("total_cases_lag0", "total_cases_lag1", "total_cases_horizon1"),
            stringsAsFactors = FALSE),
        kernel_fn = pdtmvn_kernel,
        rkernel_fn = rpdtmvn_kernel,
        theta_fixed = list(
            parameterization = "bw-diagonalized-est-eigenvalues",
            continuous_vars = c("total_cases_lag0", "total_cases_lag1", "total_cases_horizon1"),
            discrete_vars = NULL,
            discrete_var_range_fns = NULL,
            lower = c(total_cases_lag0 = -Inf, total_cases_lag1 = -Inf, total_cases_horizon1 = -Inf),
            upper = c(total_cases_lag0 = Inf, total_cases_lag1 = Inf, total_cases_horizon1 = Inf)
        ),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ))


kcde_control <- create_kcde_control(X_names = "time_index",
    y_names = "total_cases",
    time_name = "time",
    prediction_horizons = 1L,
    kernel_components = kernel_components,
    crossval_buffer = ymd("2010-01-01") - ymd("2009-01-01"),
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_type = "distribution",
    loss_args = NULL)

#tmp_file <- tempfile()
#Rprof(tmp_file, gc.profiling = TRUE, line.profiling = TRUE)

## set up parallelization
registerDoMC(cores=3)


## estimate parameters using only data up through 2014
flu_kcde_fit_orig_scale <- kcde(data = ili_national[ili_national$year <= 2014, ],
    kcde_control = kcde_control)

#Rprof(NULL)
#pd <- readProfileData(tmp_file)
#options(width = 300)
#hotPaths(pd, maxdepth = 27)


## sample from predictive distribution for first week in 2015
predictive_sample <- kcde_predict(kcde_fit = flu_kcde_fit_orig_scale,
    prediction_data = ili_national[ili_national$year == 2014 & ili_national$week %in% c(52, 53), , drop = FALSE],
    leading_rows_to_drop = 1L,
    trailing_rows_to_drop = 1L,
    additional_training_rows_to_drop = NULL,
    prediction_type = "sample",
    n = 10000L)

ggplot() +
    geom_density(aes(x = predictive_sample)) +
    geom_vline(aes(xintercept = ili_national$total_cases[ili_national$year == 2015 & ili_national$week == 1]),
        colour = "red") +
    xlab("Total Cases") +
    ylab("Predictive Density") +
    ggtitle("Realized total cases vs. one week ahead predictive density\nWeek 1 of 2015") +
    theme_bw()



## obtain kernel weights and centers
predictive_kernel_weights_and_centers <- kcde_predict(kcde_fit = flu_kcde_fit_orig_scale,
    prediction_data = ili_national[ili_national$year == 2014 & ili_national$week %in% c(52, 53), , drop = FALSE],
    leading_rows_to_drop = 1L,
    trailing_rows_to_drop = 1L,
    additional_training_rows_to_drop = NULL,
    prediction_type = "centers-and-weights")




ggplot() +
    geom_point(aes(x = predictive_kernel_weights_and_centers$centers[, 1],
            y = predictive_kernel_weights_and_centers$weights)) +
    theme_bw()


matching_predictive_value <- compute_offset_obs_vecs(data = flu_kcde_fit_orig_scale$train_data,
    vars_and_offsets = flu_kcde_fit_orig_scale$vars_and_offsets,
    time_name = flu_kcde_fit_orig_scale$kcde_control$time_name,
    leading_rows_to_drop = 1L,
    trailing_rows_to_drop = 1L,
    additional_rows_to_drop = NULL,
    na.action = flu_kcde_fit_orig_scale$kcde_control$na.action)

ggplot() +
    geom_point(aes(y = predictive_kernel_weights_and_centers$centers[, 1],
            x = matching_predictive_value$total_cases_lag0,
            alpha = predictive_kernel_weights_and_centers$weights)) +
    theme_bw()
