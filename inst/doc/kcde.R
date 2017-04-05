## ----InitialBlock, echo = FALSE------------------------------------------
library(knitr)
opts_knit$set(concordance=TRUE)

## ----FluDataLoadData, echo = TRUE----------------------------------------
library(cdcfluview)
library(dplyr)
library(lubridate)
library(ggplot2)
library(grid)
library(kcde)

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

## ----FluDataInitialPlotTotalCases, echo = TRUE---------------------------
ggplot() +
    geom_line(aes(x = as.Date(time), y = total_cases), data =
ili_national) +
    geom_vline(aes(xintercept = as.numeric(as.Date(time))),
        colour = "grey",
        data = ili_national[is.na(ili_national$total_cases), ]) +
    scale_x_date() +
    xlab("Time") +
    ylab("Total Cases") +
    theme_bw()

## ----FluDataHistogramPlotTotalCases, echo = TRUE-------------------------
hist_df <- rbind(
	data.frame(value = ili_national$total_cases,
    	variable = "Total Cases"),
    data.frame(value = log(ili_national$total_cases),
    	variable = "log(Total Cases)")
)

ggplot(aes(x = value), data = hist_df) +
    geom_histogram() +
    facet_wrap( ~ variable, ncol = 2) +
    xlab("Total Cases") +
    theme_bw()

## ----FluDataACFPlotTotalCases, echo = TRUE-------------------------------
last_na_ind <- max(which(is.na(ili_national$total_cases)))
non_na_inds <- seq(from = last_na_ind + 1, to=nrow(ili_national))
acf(ili_national$total_cases[non_na_inds],
  lag.max = 52 * 4)

## ----FluDataKernelComponentsSetup, echo = TRUE---------------------------
## Definitions of kernel components.  A couple of notes:
##   1) In the current implementation, it is required that separate kernel
##      components be used for lagged (predictive) variables and for leading
##      (prediction target) variables.
##   2) The current syntax is verbose; in a future version of the package,
##      convenience functions may be provided.

## Define kernel components -- 3 pieces:
##   1) Periodic kernel acting on time index
##   2) pdtmvn kernel acting on lagged total cases (predictive) -- all continuous
##   3) pdtmvn kernel acting on lead total cases (prediction target) -- all continuous
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
            lower = -Inf,
            upper = Inf
        ),
        theta_est = list("bw"),
        initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
        initialize_kernel_params_args = NULL,
        vectorize_kernel_params_fn = vectorize_params_pdtmvn_kernel,
        vectorize_kernel_params_args = NULL,
        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
        update_theta_from_vectorized_theta_est_args = NULL
    ))
#,
#    list(
#        vars_and_lags = vars_and_lags[3:5, ],
#        kernel_fn = pdtmvn_kernel,
#        rkernel_fn = rpdtmvn_kernel,
#        theta_fixed = NULL,
#        theta_est = list("bw"),
#        initialize_kernel_params_fn = initialize_params_pdtmvn_kernel,
#        initialize_kernel_params_args = list(
#            continuous_vars = vars_and_lags$combined_name[3:4],
#            discrete_vars = vars_and_lags$combined_name[5],
#            discrete_var_range_fns = list(
#                c_lag2 = list(a = pdtmvn::floor_x_minus_1, b = floor, in_range = pdtmvn::equals_integer, discretizer = round_up_.5))
#        ),
#        vectorize_theta_est_fn = vectorize_params_pdtmvn_kernel,
#        vectorize_theta_est_args = NULL,
#        update_theta_from_vectorized_theta_est_fn = update_theta_from_vectorized_theta_est_pdtmvn_kernel,
#        update_theta_from_vectorized_theta_est_args = list(
#            parameterization = "bw-diagonalized-est-eigenvalues"
#        )
#    ))

## ----FluDataKCDESetup, echo=TRUE-----------------------------------------
kcde_control <- create_kcde_control(X_names = "time_index",
    y_names = "total_cases",
    time_name = "time",
    prediction_horizons = 1L,
    kernel_components = kernel_components,
    crossval_buffer = ymd("2010-01-01") - ymd("2009-01-01"),
    loss_fn = neg_log_score_loss,
    loss_fn_prediction_args = list(
        prediction_type = "distribution",
        log = TRUE),
    filter_control <- NULL,
    loss_args = NULL,
    prediction_inds_not_included = c())

## ----FluDataKCDEEstimation, echo=TRUE------------------------------------
# flu_kcde_fit_orig_scale <- kcde(data = ili_national[ili_national$year <= 2014, ],
#    kcde_control = kcde_control)

## ----FluDataKCDEPredictiveSampleAndPlot, echo=TRUE-----------------------
# predictive_sample <- kcde_predict(kcde_fit = flu_kcde_fit_orig_scale,
#         prediction_data = ili_national[ili_national$year == 2014 & ili_national$week == 53, , drop = FALSE],
#         leading_rows_to_drop = 0,
#         trailing_rows_to_drop = 1L,
#         additional_training_rows_to_drop = NULL,
#         prediction_type = "sample",
#         n = 10000L)

# ggplot() +
# 	geom_density(aes(x = predictive_sample)) +
# 	geom_vline(aes(xintercept = ili_national$total_cases[ili_national$year == 2015 & ili_national$week == 1]),
# 		colour = "red") +
# 	xlab("Total Cases") +
# 	ylab("Predictive Density") +
# 	ggtitle("Realized total cases vs. one week ahead predictive density\nWeek 1 of 2015") +
# 	theme_bw()

