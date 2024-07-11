###############################################################################
#
#    sgs: Sparse-group SLOPE (Sparse-group Sorted L1 Penalized Estimation)
#    Copyright (C) 2023 Fabio Feser
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

# File contains functions for fitting a single SGS/gSLOPE model or a pathwise solution

fit_path <- function(X,y,groups, groupIDs, type, lambda_path, alpha, intercept, pen_slope_org, pen_gslope_org, X_scale, X_center, y_mean, wt, wt_per_grp, model_type,
                        num_vars, num_groups, path_length, num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose){
  machine_tol = .Machine$double.eps

  # penalty checks
  if (any(pen_slope_org < 0) | any(pen_gslope_org < 0)){
    stop("penalty sequences must be positive")
  }
  out = c()
  if (intercept){out$beta = matrix(0, nrow=num_vars+1, ncol=path_length)} else{out$beta = matrix(0, nrow=num_vars, ncol=path_length)}
  out$group_effects = matrix(0, nrow = num_groups, ncol=path_length)
  out$selected_var = list()
  out$selected_grp = list()
  out$pen_slope = pen_slope_org
  out$pen_gslope = pen_gslope_org
  out$lambda = lambda_path 
  out$type = type
  out$standardise = "l2"
  out$intercept = intercept
  out$num_it = rep(0,path_length)
  out$success = rep(0,path_length)
  out$certificate = rep(0,path_length)
  out$x = matrix(0, nrow=num_vars, ncol=path_length)
  out$z = matrix(0, nrow=num_vars, ncol=path_length)
  out$u = matrix(0, nrow=num_vars, ncol=path_length)

  for (lambda_id in 1:path_length){
    pen_slope = alpha*lambda_path[lambda_id]*pen_slope_org
    pen_gslope = (1-alpha)*lambda_path[lambda_id]*pen_gslope_org
    # -------------------------------------------------------------
    # run ATOS
    # ------------------------------------------------------------- 
    if (lambda_id == 1 | all(abs(out$x[,lambda_id-1]) <= machine_tol) | all(abs(out$u[,lambda_id-1]) <= machine_tol) ){
      warm_x = rep(0,num_vars)
      warm_u = rep(0,num_vars)
    } else {
      warm_x = out$x[,lambda_id-1]
      warm_u = out$u[,lambda_id-1]
    }

    if (type == "logistic" & intercept){
      fit_out = run_atos_log_inter(X,y, groups, groupIDs, pen_slope, pen_gslope, x0 = warm_x, u = warm_u, wt, num_vars, num_groups,
                      num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose=FALSE)
    } else {  
      fit_out = run_atos(X,y, groups, groupIDs, pen_slope, pen_gslope, x0 = warm_x, u = warm_u, wt, num_vars, num_groups,
                      num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose=FALSE)
    }
    
    if (fit_out$success == 0){ # check for convergence
      warning("algorithm did not converge, try using more iterations")
    } 

    # -------------------------------------------------------------
    # generate output
    # ------------------------------------------------------------- 
    x = fit_out$x
    z = fit_out$z
    beta_tmp = rep(0,num_vars)
    if (max((x-z)^2) < 1e-3 & mean((x-z)^2) < 1e-3){ # if solutions are very similar, pick more stable version
      if (length(which(x!=0)) <= length(which(z!=0))){ # Picking the solution with less residual values, if this is true, x is picked
        beta_tmp = as.matrix(x)
      } else {
        beta_tmp = as.matrix(z)
      }
    } else { # if solutions aren't similar, pick x
      beta_tmp = as.matrix(x)
    }

    if (type == "logistic" & intercept){
      X = X[,-1]
      log_intercept = beta_tmp[1]
      beta_tmp = beta_tmp[-1]
      groups = groups[-1]-1
      num_vars = num_vars - 1
    }
    # scale beta depending on transformations
    beta_tmp = beta_tmp/X_scale

    if (length(beta_tmp[beta_tmp!=0]) != 0){
      if (min(abs(beta_tmp[beta_tmp!=0])) < 1e-3){ # remove small resid values if solution not stable
        threshold_x = quantile(abs(beta_tmp[which(abs(beta_tmp)>(1e-4))]))[4]*1e-3
        threshold_x = quantile(abs(beta_tmp[which(abs(beta_tmp)>=threshold_x)]))[4]*1e-2
        if (!is.na(threshold_x) & lambda_path[lambda_id]!=0) { # When lambda = 0, we don't want to remove small values, as no penalisation is occuring
          threshold_x = ifelse(threshold_x>1e-2,1e-2, threshold_x) # if threshold too big, set to 1e-2 
          beta_tmp = ifelse(abs(beta_tmp)>threshold_x,beta_tmp,0)
          beta_tmp = as.matrix(beta_tmp)
        }
      }
    }

    out$selected_var[[lambda_id]] = which(beta_tmp!=0)
    which_groups_out = which_groups(beta=beta_tmp,groups=groups)
    out$selected_grp[[lambda_id]] = which_groups_out[[1]]
    out$group_effects[,lambda_id] = as.vector(which_groups_out[[2]])

    if (intercept){ # get beta back to original scale
      if (type == "linear"){
      beta_tmp = as.matrix(c(y_mean - sum(X_center*beta_tmp),beta_tmp))
      } else {
      beta_tmp = as.matrix(c(log_intercept,beta_tmp))
      }
    } 

    if (is.null(colnames(X))){ # Add variable names to output
      if (intercept){
          rownames(beta_tmp) = c("(Intercept)", paste0("v", 1:(num_vars)))
        } else {
          rownames(beta_tmp) = paste0("v", 1:num_vars)
        }
      } else {
      if (intercept){
          rownames(beta_tmp) = c("(Intercept)", colnames(X))
        } else {
          rownames(beta_tmp) = colnames(X)
        }
    }
    out$beta[,lambda_id] = beta_tmp
    out$z[,lambda_id] = z
    out$x[,lambda_id] = x
    out$u[,lambda_id] = fit_out$u
    out$success[lambda_id] = fit_out$success
    out$num_it[lambda_id] = fit_out$it
    out$certificate[lambda_id] = fit_out$certificate

    if (verbose){print(paste0("Lambda ", lambda_id,"/",path_length, " done"))}
  }
  return(out)
}

fit_one <- function(X,y,groups, groupIDs, type, lambda, alpha, intercept, pen_slope_org, pen_gslope_org, x0=rep(0,ncol(X)), u=rep(0,ncol(X)), X_scale, X_center, 
                        y_mean=rep(0,nrow(X)), wt, wt_per_grp, model_type, num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose){
  # penalty checks
  if (any(pen_slope_org < 0) | any(pen_gslope_org < 0)){
    stop("penalty sequences must be positive")
  }

  num_groups = length(unique(groups))
  num_vars = length(groups)
  # weights
  pen_slope = alpha*lambda*pen_slope_org
  pen_gslope = (1-alpha)*lambda*pen_gslope_org

  # -------------------------------------------------------------
  # run ATOS
  # ------------------------------------------------------------- 
  if (type == "logistic" & intercept){ # doesn't work
    fit_out = run_atos_log_inter(X,y, groups, groupIDs, pen_slope, pen_gslope, x0, u, wt, num_vars, num_groups,
                      num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose)
  } else { # run SGS
    fit_out = run_atos(X,y, groups, groupIDs, pen_slope, pen_gslope, x0, u, wt, num_vars, num_groups,
                      num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose)
  }

  if (fit_out$success == 0){ # check for convergence
    warning("algorithm did not converge, try using more iterations")
  } else {if (verbose==TRUE){print("Algorithm converged")}}

  # -------------------------------------------------------------
  # generate output
  # ------------------------------------------------------------- 
  out = c()
  out$beta = rep(0,num_vars)
  x = fit_out$x
  z = fit_out$z
  if (max((x-z)^2) < 1e-3 & mean((x-z)^2) < 1e-3){ # if solutions are very similar, pick more stable version
    if (length(which(x!=0)) <= length(which(z!=0))){ # Picking the solution with less residual values, if this is true, x is picked
      beta_tmp = as.matrix(x)
    } else {
      beta_tmp = as.matrix(z)
    }
  } else { # if solutions aren't similar, pick x
    beta_tmp = as.matrix(x)
  }

  if (type == "logistic" & intercept){
    X = X[,-1]
    log_intercept = beta_tmp[1]
    beta_tmp = beta_tmp[-1]
    groups = groups[-1]-1
    num_vars = num_vars - 1
  }
  # scale beta depending on transformations
  beta_tmp = beta_tmp/X_scale

  if (length(beta_tmp[beta_tmp!=0]) != 0){
    if (min(abs(beta_tmp[beta_tmp!=0])) < 1e-3){ # remove small resid values if solution not stable
      threshold_x = quantile(abs(beta_tmp[which(abs(beta_tmp)>(1e-4))]))[4]*1e-3
      threshold_x = quantile(abs(beta_tmp[which(abs(beta_tmp)>=threshold_x)]))[4]*1e-2
      if (!is.na(threshold_x) & lambda!=0) { # When lambda = 0, we don't want to remove small values, as no penalisation is occuring
        threshold_x = ifelse(threshold_x>1e-2,1e-2, threshold_x) # if threshold too big, set to 1e-2 
        beta_tmp = ifelse(abs(beta_tmp)>threshold_x,beta_tmp,0)
        beta_tmp = as.matrix(beta_tmp)
      }
    }
  }
  which_groups_out = which_groups(beta=beta_tmp,groups=groups)
  out$group_effects = which_groups_out[[2]]
  out$selected_var = which(beta_tmp!=0)
  out$selected_grp = which_groups_out[[1]]

  if (intercept){ # get beta back to original scale
    if (type == "linear"){
      beta_tmp = as.matrix(c(y_mean - sum(X_center*beta_tmp),beta_tmp))
    } else {
      beta_tmp = as.matrix(c(log_intercept,beta_tmp))
    }
  } 

  if (is.null(colnames(X))){ # Add variable names to output
    if (intercept){
        rownames(beta_tmp) = c("(Intercept)", paste0("v", 1:(num_vars)))
      } else {
        rownames(beta_tmp) = paste0("v", 1:num_vars)
      }
    } else {
    if (intercept){
        rownames(beta_tmp) = c("(Intercept)", colnames(X))
      } else {
        rownames(beta_tmp) = colnames(X)
      }
  }
  out$beta = beta_tmp
  out$pen_slope = pen_slope_org
  out$pen_gslope = pen_gslope_org
  out$lambda=lambda
  out$type = type
  out$standardise = "l2"
  out$intercept = intercept
  out$num_it = fit_out$it
  out$success = fit_out$success
  out$certificate = fit_out$certificate
  out$x = x
  out$z = z
  out$u = fit_out$u
  return(out)
}