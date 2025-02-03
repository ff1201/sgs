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

# File contains ATOS fitting code applied specifically to SGS/gSLOPE

run_atos <- function(X, y, groups, groupIDs, pen_slope, pen_gslope, x0, u, wt, num_vars, num_groups,
                    num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose){
  # set values
  inv_wt = 1/wt
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  LS_EPS = .Machine$double.eps # R accuracy

  # initial fitting values
  tX = Matrix::t(X)
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, mult_fcn=mult_fcn, x0=x0, X=X,y=y,num_obs=num_obs, tX=tX, crossprod_mat=crossprod_mat)
  z = proxGroupSortedL1(y=x0*wt, lambda=pen_gslope*step_size,group=groups,group_id=groupIDs,num_groups=num_groups)
  z = inv_wt*z
  Xbeta = mult_fcn(X,z)
  fz = f(y, Xbeta, num_obs, crossprod_mat)
  grad_fz = mult_fcn(tX,f_grad(y, Xbeta, num_obs)) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  x = sortedL1Prox(x=z - (step_size * (grad_fz)) ,lambda=pen_slope*step_size)

  # fitting
  for (it in 1:max_iter){
    Xbeta = mult_fcn(X,z)
    fz = f(y, Xbeta, num_obs, crossprod_mat)
    grad_fz = mult_fcn(tX,f_grad(y, Xbeta, num_obs))
    x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))), lambda=pen_slope*step_size, method="stack")
    incr = x - z
    norm_incr = norm_vec(incr)
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))), lambda=pen_slope*step_size, method="stack")
        incr = x - z
        norm_incr = norm_vec(incr)
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y, mult_fcn(X,x), num_obs, crossprod_mat) - rhs        
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     

    z = proxGroupSortedL1(y=wt*x + (step_size/wt)*u, lambda=pen_gslope*step_size,group=groups,group_id=groupIDs,num_groups=num_groups)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size
    
    if (certificate < tol){ # Check for convergence
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = ifelse(it==max_iter,0,1)
out$certificate = certificate
out$it = it
return(out)
}

run_atos_log_inter <- function(X, y, groups, groupIDs, pen_slope, pen_gslope, x0, u, wt, num_vars, num_groups,
                    num_obs, max_iter, backtracking, max_iter_backtracking, f, f_grad, mult_fcn, crossprod_mat, tol, verbose){ # needs updating
  # set values
  inv_wt = 1/wt
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  LS_EPS = .Machine$double.eps # R accuracy
  tX = Matrix::t(X)
  # initial fitting values
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, mult_fcn=mult_fcn, x0=x0, X=X,y=y,num_obs=num_obs, tX=tX, crossprod_mat=crossprod_mat)
  prox_input = x0*wt
  z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups,group_id=groupIDs, num_groups)
  z = c(prox_input[1],z)
  z = inv_wt*z
  Xbeta = mult_fcn(X,z)
  fz = f(y, Xbeta, num_obs, crossprod_mat)
  grad_fz = f_grad(y, Xbeta, num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  prox_input = z - (step_size * (grad_fz)) 
  x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
  x = c(prox_input[1],x)
  # fitting
  for (it in 1:max_iter){
    Xbeta = mult_fcn(X,z)
    fz = f(y, Xbeta, num_obs, crossprod_mat)
    grad_fz = f_grad(y, Xbeta, num_obs) # loss gradient at z
    prox_input = z - (step_size * (u + (grad_fz))) 
    x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
    x = c(prox_input[1],x)
    incr = x - z
    norm_incr = norm_vec(incr)
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        prox_input = z - (step_size * (u + (grad_fz)))
        x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
        x = c(prox_input[1],x)
        incr = x - z
        norm_incr = norm_vec(incr)
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y, mult_fcn(X,x), num_obs, crossprod_mat) - rhs  
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     
    prox_input = wt*x + (step_size/wt)*u
    z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups,group_id=groupIDs, num_groups=num_groups)
    z = c(prox_input[1],z)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size

    if (certificate < tol){ # Check for convergence
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = ifelse(it==max_iter,0,1)
out$certificate = certificate
out$it = it
return(out)
}
