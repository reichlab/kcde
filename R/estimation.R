## Functions for parameter estimation
##
## kcde
## est_kcde_params_stepwise_crossval
## est_kcde_params_stepwise_crossval_one_potential_step
## kcde_crossval_estimate_parameter_loss

#' Estimate the parameters for kcde.  There is redundancy here in that X_names,
#' y_names, time_name are all included in the kcde_control object as well as
#' parameters to this function.  Decide on an interface.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' @param kcde_control a list of parameters kcde_controlling how the fit is done.
#'     See the documentation for kcde_control.
#' 
#' @return an object representing an estimated kcde model.  Currently, a list
#'     with 7 components.
kcde <- function(X_names,
        y_names,
        time_name,
        data,
        kcde_control) {
    ## get/validate kcde_control argument
    if(missing(kcde_control)) {
        kcde_control <- create_kcde_control_default(X_names, y_names, time_name, data)
        warning("kcde_control argument not supplied to kcde -- using defaults, which may be bad")
    } else {
        validate_kcde_control(kcde_control, X_names, y_names, time_name, data)
    }
    
    ## estimate lags and kernel parameters via cross-validation
    param_estimates <- est_kcde_params_stepwise_crossval(data, kcde_control)
    
    return(list(kcde_control=kcde_control,
        X_names=X_names,
        y_names=y_names,
        time_name=time_name,
		vars_and_lags=param_estimates$vars_and_lags,
        theta_hat=param_estimates$theta_hat,
        train_data=data))
}

#' Use a stepwise procedure to estimate parameters and select lags by optimizing
#' a cross-validation estimate of predictive performance
#' 
#' @param data the data frame to use in performing cross validation
#' @param kcde_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with two components: vars_and_lags is the estimated "optimal" lags
#'     to use for each variable, and theta_hat is the estimated "optimal"
#'     kernel parameters to use for each combination of variable and lag
est_kcde_params_stepwise_crossval <- function(data, kcde_control) {
    all_vars_and_lags <- rbind.fill(lapply(kcde_control$kernel_components,
		function(kernel_component) {
			kernel_component$vars_and_lags
		}))
    
    ## initialize cross-validation process
    ## the variable selected in previous iteration
    selected_var_lag_ind <- NULL
    ## cross-validation estimate of loss associated with current estimates
    crossval_prediction_loss <- Inf
    ## initial parameters: no variables/lags selected, no kernel parameters
    vars_and_lags <- data.frame(var_name = NA,
		lag_value = NA,
		combined_name = NA)
    theta_hat <- vector("list", length(kcde_control$kernel_components))
    
    all_evaluated_models <- list()
    all_evaluated_model_descriptors <- list()
    
    repeat {
        ## get cross-validation estimates of performance for model obtained by
        ## adding or removing each variable/lag combination (except for the one
        ## updated in the previous iteration) from the model, as well as the
        ## corresponding parameter estimates
        
        ## commented out use of foreach for debugging purposes
#        crossval_results <- foreach(i=seq_len(nrow(all_vars_and_lags)),
#            .packages=c("kcde", kcde_control$par_packages),
#            .combine="c") %dopar% {
        crossval_results <- lapply(seq_len(nrow(all_vars_and_lags)),
			function(i) {
	        	descriptor_current_model <- update_vars_and_lags(vars_and_lags,
					all_vars_and_lags[i, "var_name"],
					all_vars_and_lags[i, "lag_value"])
        		
	        	model_i_previously_evaluated <- any(sapply(
					all_evaluated_model_descriptors,
	        		function(descriptor) {
	        			identical(descriptor_current_model, descriptor)
	        		}
				))
        		
	            if(!model_i_previously_evaluated) {
	            	potential_step_result <- 
	                    est_kcde_params_stepwise_crossval_one_potential_step(
	                        prev_vars_and_lags=vars_and_lags,
	                        prev_theta=theta_hat,
	                        update_var_name=all_vars_and_lags[i, "var_name"],
	                        update_lag_value=all_vars_and_lags[i, "lag_value"],
	                        data=data,
	                        kcde_control=kcde_control)
                
	                all_evaluated_models <-
	                	c(all_evaluated_models,
	                		list(potential_step_result))
	               	all_evaluated_model_descriptors <-
	               		c(all_evaluated_model_descriptors,
	               			potential_step_result["all_vars_and_lags"])
	                
	                return(potential_step_result)
#	               return(list(potential_step_result)) # put results in a list so that combine="c" is useful
	            } else {
	            	return(NULL)
	            }
#	        }
	        }
		)
       
        ## drop elements corresponding to previously explored models
        non_null_components <- sapply(crossval_results,
        	function(component) { !is.null(component) }
        )
        crossval_results <- crossval_results[non_null_components]
        
        ## pull out loss achieved by each model, find the best value
        loss_achieved <- sapply(crossval_results, function(component) {
            component$loss
        })
        optimal_loss_ind <- which.min(loss_achieved)
       
#        print("loss achieved is:")
#        print(loss_achieved)
#        print("\n")
        
        ## either update the model and keep going or stop the search
        if(loss_achieved[optimal_loss_ind] < crossval_prediction_loss) {
            ## found a model improvement -- update and continue
            selected_var_lag_ind <- optimal_loss_ind
            crossval_prediction_loss <- loss_achieved[selected_var_lag_ind]
			vars_and_lags <-
				crossval_results[[selected_var_lag_ind]]$vars_and_lags
            theta_hat <- crossval_results[[selected_var_lag_ind]]$theta
        } else {
            ## could not find a model improvement -- stop search
            break
        }
    }

    return(list(vars_and_lags=vars_and_lags,
        theta_hat=theta_hat))
}


#' Estimate the parameters theta and corresponding cross-validated estimate of
#' loss for one possible model obtained by adding or removing a variable/lag
#' combination from the model obtained in the previous iteration of the stepwise
#' search procedure.
#' 
#' @param prev_vars_and_lags list representing combinations of variables and lags
#'     included in the model obtained at the previous step
#' @param prev_theta list representing the kernel parameter estimates obtained
#'     at the previous step
#' @param update_var_name the name of the variable to try adding or removing
#'     from the model
#' @param update_lag_value the value of the lag for the variable specified by
#'     update_var_name to try adding or removing from the model
#' @param data the data frame with observations used in estimating model
#'     parameters
#' @param kcde_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with three components: loss is a cross-validation estimate of
#'     the loss associated with the estimated parameter values for the given
#'     model, lags is a list representing combinations of variables and lags
#'     included in the updated model, and theta is a list representing the
#'     kernel parameter estimates in the updated model
est_kcde_params_stepwise_crossval_one_potential_step <- function(
		prev_vars_and_lags,
    	prev_theta,
    	update_var_name,
    	update_lag_value,
    	data,
    	kcde_control) {
    ## updated variable and lag combinations included in model
    updated_vars_and_lags <- update_vars_and_lags(prev_vars_and_lags,
		update_var_name,
		update_lag_value)
    
    ## initial values for theta; a list with two components:
    ##     theta_est, vector of initial values in vector form on scale appropriate for
    ##         estimation.
    ##     theta_fixed, list of values that will not be estimated, one component
    ##         for each component kernel function
    theta_init <- initialize_theta(prev_theta,
    	update_var_name,
    	update_lag_value,
    	data,
    	kcde_control)
    	
    theta_est_init <- extract_vectorized_theta_est_from_theta(theta_init,
    	kcde_control)
    
    ## optimize parameter values
    optim_result <- optim(par=theta_est_init,
        fn=kcde_crossval_estimate_parameter_loss,
#        gr = gradient_kcde_crossval_estimate_parameter_loss,
		gr=NULL,
		theta=theta_init,
        vars_and_lags=updated_vars_and_lags,
        data=data,
        kcde_control=kcde_control,
        method="L-BFGS-B",
        lower=-50,
        #		upper=10000,
        #control=list(),
        hessian=FALSE)
    
#    print(optim_result)
    
    ## convert back to list and original parameter scale
    updated_theta <- update_theta_from_vectorized_theta_est(optim_result$par,
        theta_init,
        kcde_control)
    
    return(list(
        loss=optim_result$value,
        vars_and_lags=updated_vars_and_lags,
        theta=updated_theta
    ))
}

#' Using cross-validation, estimate the loss associated with a particular set
#' of lags and kernel function parameters.
#' 
#' @param theta_est_vector vector of kernel function parameters that are being
#'     estimated
#' @param theta list of kernel function parameters, both those that are being
#'     estimated and those that are out of date.  Possibly the values of
#'     parameters being estimated are out of date; they will be replaced with
#'     the values in theta_est_vector.
#' @param vars_and_lags list representing combinations of variables and lags
#'     included in the model
#' @param data the data frame to use in performing cross validation
#' @param kcde_control a list of parameters specifying how the fitting is done
#' 
#' @return numeric -- cross-validation estimate of loss associated with the
#'     specified parameters
kcde_crossval_estimate_parameter_loss <- function(theta_est_vector,
		theta,
		vars_and_lags,
    	data,
    	kcde_control) {
	max_lag <- max(unlist(lapply(kcde_control$kernel_components,
		function(kernel_component) {
			kernel_component$vars_and_lags$lag_value
		}
	)))
    ## set up theta list containing both the kernel parameters that are being
    ## estimated and the kernel parameters that are being held fixed
    ## also, transform back to original scale
    ## convert back to list and original parameter scale
    theta <- update_theta_from_vectorized_theta_est(theta_est_vector,
        theta,
        kcde_control)
    
    ## create data frame of "examples" -- lagged observation vectors and
    ## corresponding prediction targets
    cross_validation_examples <- assemble_training_examples(data,
		vars_and_lags,
        kcde_control$y_names,
        leading_rows_to_drop=max_lag,
        additional_rows_to_drop=NULL,
        prediction_horizons=kcde_control$prediction_horizons,
        drop_trailing_rows=TRUE)
    
    ## This could be made more computationally efficient by computing
    ## kernel values for all relevant combinations of lags for each variable,
    ## then combining as appropriate -- currently, the same kernel value may be
    ## computed multiple times in the call to kcde_predict_given_lagged_obs
    crossval_loss_by_time_ind <- sapply(
        seq_len(nrow(cross_validation_examples$lagged_obs)),
        function(t_pred) {
            ## get training indices -- those indices not within
            ## t_pred +/- kcde_control$crossval_buffer
            t_train <- seq_len(nrow(cross_validation_examples$lagged_obs))
            t_train <- t_train[!(t_train %in%
                seq(from=t_pred - kcde_control$crossval_buffer,
                    to=t_pred + kcde_control$crossval_buffer))]
            
            ## calculate kernel weights and centers for prediction at
            ## prediction_lagged_obs based on train_lagged_obs and
            ## train_lead_obs
            ## assemble lagged and lead observations -- subsets of
            ## cross_validation_examples given by t_pred and t_train
            ## we can re-use the weights at different prediction_target_inds,
            ## and just have to adjust the kernel centers
            prediction_target_ind <- 1
            
            train_lagged_obs <- cross_validation_examples$lagged_obs[
                t_train, , drop=FALSE]
            train_lead_obs <- cross_validation_examples$lead_obs[
                t_train, prediction_target_ind, drop=FALSE]
            prediction_lagged_obs <- 
                cross_validation_examples$lagged_obs[
                    t_pred, , drop=FALSE]
            prediction_lead_obs <-
                cross_validation_examples$lead_obs[
                    t_pred, prediction_target_ind, drop=FALSE]
            
            ## for each prediction target variable, compute loss
            crossval_loss_by_prediction_target <- sapply(
                seq_len(ncol(cross_validation_examples$lead_obs)),
                function(prediction_target_ind) {
                    ## calculate and return value of loss function based on
                    ## prediction and realized value
                    loss_args <- kcde_control$loss_fn_args
                    loss_args$prediction_result <- kcde_predict_given_lagged_obs(
		                train_lagged_obs=train_lagged_obs,
		                train_lead_obs=train_lead_obs,
		                prediction_lagged_obs=prediction_lagged_obs,
		                prediction_lead_obs=prediction_lead_obs,
		                kcde_fit=list(theta_hat=theta,
		                    kcde_control=kcde_control
		                ),
		                prediction_type=kcde_control$loss_fn_prediction_type)
                    
                    loss_args$obs <- as.numeric(
                        cross_validation_examples$lead_obs[
                            t_pred, prediction_target_ind]
                    )
                    
                    return(do.call(kcde_control$loss_fn, loss_args))
                })
            
            return(sum(crossval_loss_by_prediction_target))
        })
    
    if(any(is.na(crossval_loss_by_time_ind))) {
        ## parameters resulted in numerical instability
        stop("NAs in cross validated estimate of loss")
        ## old solution was to return largest non-infinite value
        return(.Machine$double.xmax)
    } else {
        return(sum(crossval_loss_by_time_ind))
    }
}
