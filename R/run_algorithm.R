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

run_atos_log_inter <- function(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, x0, u,
                          pen_slope, pen_gslope, f, f_grad,f_opts, f_grad_opts, crossprod_mat,verbose,groups_pen){
  # set values
  wt = as.numeric(sqrt(rep(table(groups),table(groups)))) # Adjust for group weights
  inv_wt = 1/wt
  group_ids = getGroupID(groups) 
  group_ids_pen = getGroupID(groups_pen)
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  num_groups = length(unique(groups))
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  # initial fitting values
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
  prox_input = x0*wt
  z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups_pen,group_id=group_ids_pen)
  z = c(prox_input[1],z)
  z = inv_wt*z
  fz = f(y=y, X=X, input = z, num_obs = num_obs)
  grad_fz = f_grad(y=y, X=X, input = z, num_obs = num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  prox_input = z - (step_size * (grad_fz)) 
  x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
  x = c(prox_input[1],x)
  # fitting
  for (it in 1:max_iter){
    fz = f(y=y, X=X, input = z, num_obs = num_obs)
    grad_fz = f_grad(y=y, X=X, input =z, num_obs = num_obs)
    prox_input = z - (step_size * (u + (grad_fz))) 
    x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
    x = c(prox_input[1],x)
    incr = x - z
    norm_incr = norm(incr,type="2")
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        prox_input = z - (step_size * (u + (grad_fz)))
        x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
        x = c(prox_input[1],x)
        incr = x - z
        norm_incr = norm(incr,type="2")
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y=y, X=X, input = x, num_obs = num_obs) - rhs        
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     
    prox_input = wt*x + (step_size/wt)*u
    z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups_pen,group_id=group_ids_pen)
    z = c(prox_input[1],z)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size

    if (certificate < tol){ # Check for convergence
      success = 1
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = success
out$certificate = certificate
out$it = it
return(out)
}

run_atos <- function(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, x0, u,
                          pen_slope, pen_gslope, f, f_grad, f_opts, f_grad_opts, crossprod_mat,verbose){
  # set values
  wt = as.numeric(sqrt(rep(table(groups),table(groups)))) # Adjust for group weights
  inv_wt = 1/wt
  group_ids = getGroupID(groups) 
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  num_groups = length(unique(groups))
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  # initial fitting values
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
  z = proxGroupSortedL1(y=x0*wt, lambda=pen_gslope*step_size,group=groups,group_id=group_ids)
  z = inv_wt*z
  fz = f(y=y, X=X, input = z, num_obs = num_obs)
  grad_fz = f_grad(y=y, X=X, input = z, num_obs = num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  x = sortedL1Prox(x=z - (step_size * (grad_fz)) ,lambda=pen_slope*step_size)

  # fitting
  for (it in 1:max_iter){
    fz = f(y=y, X=X, input = z, num_obs = num_obs)
    grad_fz = f_grad(y=y, X=X, input =z, num_obs = num_obs)
    x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))) ,lambda=pen_slope*step_size)
    incr = x - z
    norm_incr = norm(incr,type="2")
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))) ,lambda=pen_slope*step_size)
        incr = x - z
        norm_incr = norm(incr,type="2")
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y=y, X=X, input = x, num_obs = num_obs) - rhs        
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     

    z = proxGroupSortedL1(y=wt*x + (step_size/wt)*u, lambda=pen_gslope*step_size,group=groups,group_id=group_ids)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size

    if (certificate < tol){ # Check for convergence
      success = 1
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = success
out$certificate = certificate
out$it = it
return(out)
}