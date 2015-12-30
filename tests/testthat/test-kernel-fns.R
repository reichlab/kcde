library(kcde)
library(magrittr)
library(plyr)
library(dplyr)

## all functions in kernel-fns.R
##   leading P means test written and passing
##   leading F means test written but failing
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
##
## S periodic_kernel
##   get_col_inds_continuous_dicrete_vars_used
##   pdtmvn_kernel
## S compute_pdtmvn_kernel_bw_params_from_bw_eigen
##   vectorize_params_pdtmvn_kernel
##   update_theta_from_vectorized_theta_est_pdtmvn_kernel
## S initialize_params_pdtmvn_kernel

context("kernel-fns -- no tests implemented")

#test_that("initialize_params_pdtmvn_kernel works", {
#	parameterization <- "bw-diagonalized-est-eigenvalues"
#
#	n <- 100
#	x <- matrix(rnorm(3*n), nrow = n)
#	x[c(2, 10:95), 3] <- floor(x[c(3, 10:95), 1])
#	colnames(x) <- letters[24:26]
#	
#	continuous_vars <- letters[24:25]
#	discrete_vars <- letters[26]
#	
#	kcde_control <- list(junk = 1:5)
#	
#	actual <- initialize_params_pdtmvn_kernel(parameterization,
#		x,
#		continuous_vars,
#		discrete_vars,
#		kcde_control)
#	
#	temp <- covRob(x)$cov
#})
#
#test_that("pdtmvn_kernel")