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


last_na_ind <- max(which(is.na(ili_national$total_cases)))
ili_national <- ili_national[seq(from = last_na_ind + 1, to = nrow(ili_national)), , drop=FALSE]

temp <- acf(ili_national$total_cases, type = "covariance", plot = FALSE, lag.max = 3)
temp <- temp$acf[, , 1]

orig_cov <- matrix(NA, nrow = 4, ncol = 4)
junk <- apply(expand.grid(1:4, 1:4), 1, function(inds) {
	ind1 <- inds[1]
	ind2 <- inds[2]
	lag <- abs(ind2 - ind1)
	orig_cov[ind1, ind2] <<- temp[lag + 1]
})



k <- 10
smooth_total_cases <- stats::filter(ili_national$total_cases, rep(1/k, k), sides = 1)

temp_smooth <- acf(smooth_total_cases[!is.na(smooth_total_cases)], type = "covariance", plot = FALSE, lag.max = 3)
temp_smooth <- temp_smooth$acf[, , 1]

smooth_cov <- matrix(NA, nrow = 4, ncol = 4)
junk <- apply(expand.grid(1:4, 1:4), 1, function(inds) {
	ind1 <- inds[1]
	ind2 <- inds[2]
	lag <- abs(ind2 - ind1)
	smooth_cov[ind1, ind2] <<- temp_smooth[lag + 1]
})

eigen(smooth_cov)
eigen(orig_cov)

