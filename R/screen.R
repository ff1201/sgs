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

# File contains screening code 

screen_strong <- function(X, y, groups, groupIDs, type, lambda_path, alpha, pen_slope_org, pen_gslope_org, X_scale, num_vars, wt, path_length, model, 
                          var_screen_fcn, grp_screen_fcn, kkt_fcn, num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose){
  # -------------------------------------------------------------
  # initial set-up
  # ------------------------------------------------------------- 
  screen_fitting_options = list(num_obs=num_obs, max_iter=max_iter, backtracking=backtracking, max_iter_backtracking = max_iter_backtracking,
                          f=f,f_grad=f_grad,mult_fcn=mult_fcn,crossprod_mat=crossprod_mat,tol=tol,verbose=FALSE)
                          
  machine_tol = .Machine$double.eps
  out = c()
  out$screen_set_var = list()
  out$screen_set_var[[1]] = "no screening performed"
  out$screen_set_grp = list()
  out$screen_set_grp[[1]] = "no screening performed"
  out$kkt_violations_var = list()
  out$kkt_violations_var[[1]] = "no screening performed"
  out$kkt_violations_grp = list()
  out$kkt_violations_grp[[1]] = "no screening performed"
  out$epsilon_set_var = list()
  out$epsilon_set_var[[1]] = "no screening performed"
  out$epsilon_set_grp = list()
  out$epsilon_set_grp[[1]] = "no screening performed"
  out$active_set_var = list()
  out$active_set_grp = list()
  out$num_it = rep(0,path_length)
  out$success = rep(0,path_length)
  out$certificate = rep(0,path_length)
  out$x = matrix(0,nrow=num_vars,ncol=path_length)
  out$z = matrix(0,nrow=num_vars,ncol=path_length)
  out$u = matrix(0,nrow=num_vars,ncol=path_length)
  tbl_grps = table(groups)
  tbl_grps_sqrt = sqrt(table(groups))
  tX = Matrix::t(X)

  # -------------------------------------------------------------
  # Fit model for lambda_max
  # ------------------------------------------------------------- 
  warm_x0 = rep(0,num_vars)
  warm_u = rep(0,num_vars)
  out$beta = matrix(0,nrow=num_vars,ncol=path_length)
  out$group_effects = matrix(0,nrow=length(tbl_grps),ncol=path_length)
  current_model = do.call(fit_one, c(list(X,y,groups, groupIDs, type, lambda_path[1], alpha=alpha, FALSE, pen_slope_org, pen_gslope_org, x0 = warm_x0, u = warm_u, X_scale = X_scale, X_center=rep(0,num_vars),
  y_mean=rep(0,num_obs), wt=wt, wt_per_grp=tbl_grps_sqrt, model), screen_fitting_options))
  out$beta[,1] = as.vector(current_model$beta)
  current_beta = out$beta[,1]*X_scale
  out$active_set_grp[[1]] = current_model$selected_grp
  out$active_set_var[[1]] = which(current_beta!=0)
  if (any(abs(current_model$x) > machine_tol) & any(abs(current_model$u) > machine_tol)){
      warm_x0 = current_model$x
      warm_u = current_model$u     
  }

  if (verbose){print(paste0("Lambda ", 1,"/",path_length, " done"))}
   
  # -------------------------------------------------------------
  # begin screening
  # ------------------------------------------------------------- 
  for (lambda_id in 2:path_length){
    # calculate gradient
    active_set_var = out$active_set_var[[lambda_id-1]]
    active_set_grp = out$active_set_grp[[lambda_id-1]] 
    grad_vec = mult_fcn(tX,f_grad(y, mult_fcn(X,current_beta), num_obs))

    # group screening
    screen_set_grp = grp_screen_fcn(grad_vec, current_beta, tbl_grps, groupIDs, alpha, pen_slope_org, pen_gslope_org, lambda_path[lambda_id], lambda_path[lambda_id-1], tbl_grps_sqrt)

    # variable screening - skip if gslope
    if (length(screen_set_grp)>0 & model != "gslope"){ 
      screen_set_var = var_screen_fcn(grad_vec, groupIDs, screen_set_grp, alpha, pen_slope_org, lambda_path[lambda_id], lambda_path[lambda_id-1], active_set_var)
      epsilon_set_var = sort(union(active_set_var,screen_set_var)) # corresponds to capital epsilon
      epsilon_set_grp = unique(groups[epsilon_set_var]) # corresponds to capital epsilon
    } else if (model == "gslope"){
      screen_set_var = 1
      epsilon_set_grp = sort(union(active_set_grp,screen_set_grp)) # corresponds to capital epsilon
      epsilon_set_var = as.vector(unlist(groupIDs[epsilon_set_grp]))
    } else {
      screen_set_var = screen_set_grp
      epsilon_set_var = active_set_var
      epsilon_set_grp = active_set_grp
    }

    # initial fit
    kkt_violations_var = 1  
    out$kkt_violations_var[[lambda_id]] = 0
    out$kkt_violations_grp[[lambda_id]] = 0
    while (length(kkt_violations_var) > 0){
      current_beta = rep(0,num_vars)
      fitting_var = epsilon_set_var
      if (length(fitting_var) == 0){
        current_beta = rep(0,num_vars)
        current_model = c()
        current_model$selected_grp = numeric(0)
        current_model$selected_var = integer(0)
        current_model$num_it = 0
        current_model$success = 1
        current_model$certificate = 0
        current_model$x = rep(0,num_vars)
        current_model$u = rep(0,num_vars)
        current_model$z = rep(0,num_vars)
      } else {
        current_model = do.call(fit_one, c(list(X[,fitting_var],y,groups[fitting_var], getGroupID(groups[fitting_var]), type, lambda_path[lambda_id], alpha=alpha, FALSE, pen_slope_org[1:length(fitting_var)], 
                                        pen_gslope_org[1:length(epsilon_set_grp)], x0 = warm_x0[fitting_var], u = warm_u[fitting_var], X_scale = X_scale[fitting_var], X_center=rep(0,length(fitting_var)),
                                        y_mean=rep(0,num_obs), wt=wt[fitting_var], wt_per_grp=tbl_grps_sqrt[epsilon_set_grp], model), screen_fitting_options))
        current_beta[fitting_var] = as.vector(current_model$beta)*X_scale[fitting_var]
        if (any(abs(current_model$x) > machine_tol) & any(abs(current_model$u) > machine_tol)){
          warm_x0[fitting_var] = current_model$x
          warm_u[fitting_var] = current_model$u     
        }    
      } 

      # kkt check
      grad_vec = mult_fcn(tX,f_grad(y, mult_fcn(X,current_beta), num_obs))
      kkt_set = kkt_fcn(grad_vec,current_beta,groups,groupIDs,alpha,pen_slope_org,pen_gslope_org,lambda_path[lambda_id],tbl_grps,machine_tol,epsilon_set_var, current_model$selected_grp)

      # check for violations
      kkt_violations_grp = which(kkt_set[[2]]==1)[which(!(which(kkt_set[[2]]==1) %in% epsilon_set_grp))]
      if (model == "gslope"){
        kkt_violations_var = kkt_violations_grp
        epsilon_set_grp = sort(union(epsilon_set_grp,kkt_violations_grp))
        epsilon_set_var = as.vector(unlist(groupIDs[epsilon_set_grp]))
      } else {
        kkt_violations_var = which(kkt_set[[1]]==1)[which(!(which(kkt_set[[1]]==1) %in% epsilon_set_var))]
        epsilon_set_var = sort(union(epsilon_set_var,kkt_violations_var))
        epsilon_set_grp = unique(groups[epsilon_set_var]) # corresponds to capital epsilon
      }
      out$kkt_violations_var[[lambda_id]] = sort(union(kkt_violations_var,out$kkt_violations_var[[lambda_id]]))
      out$kkt_violations_grp[[lambda_id]] = sort(union(kkt_violations_grp,out$kkt_violations_grp[[lambda_id]]))
    }
    if (length(out$kkt_violations_var[[lambda_id]]) > 1){
      out$kkt_violations_var[[lambda_id]] = out$kkt_violations_var[[lambda_id]][-1]
    }
    if (length(out$kkt_violations_grp[[lambda_id]]) > 1){
      out$kkt_violations_grp[[lambda_id]] = out$kkt_violations_grp[[lambda_id]][-1]
    }

    # prepare output
    out$beta[fitting_var,lambda_id] = as.vector(current_model$beta)
    out$group_effects[epsilon_set_grp,lambda_id] = current_model$group_effects
    out$epsilon_set_var[[lambda_id]] = epsilon_set_var
    out$epsilon_set_grp[[lambda_id]] = epsilon_set_grp
    out$active_set_grp[[lambda_id]] = current_model$selected_grp
    out$active_set_var[[lambda_id]] = which(current_beta!=0)
    out$screen_set_grp[[lambda_id]] = sort(screen_set_grp)
    out$screen_set_var[[lambda_id]] = sort(screen_set_var)
    out$num_it[lambda_id] = current_model$num_it
    out$success[lambda_id] = current_model$success
    out$certificate[lambda_id] = current_model$certificate
    out$x[,lambda_id][fitting_var] = current_model$x
    out$z[,lambda_id][fitting_var] = current_model$z
    out$u[,lambda_id][fitting_var] = current_model$u
  
    if (verbose){print(paste0("Lambda ", lambda_id,"/",path_length, " done"))}
  }
  return(out)
}