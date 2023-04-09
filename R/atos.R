###############################################################################
#
#    sgs: Sparse-group SLOPE (Sparse-group Sorted L1 Penalized Estimation)
#    Copyright (C) 2022 Fabio Feser
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

#' Adaptive three operator splitting (ATOS)
#'
#' General function to perform ATOS with penalties other than SGS.
#'
#' Performs ATOS for other penalties. An example for the sparse-group lasso (SGL) is given. 
#' The algorithm is not symmetrical, but usually the difference between variations are only small numerical values, which are filtered out.
#' However, both variations should be checked regardless.
#'
#' @param X Input matrix of dimensions \eqn{p \times n}.
#' @param y Output vector of dimension \eqn{n}.
#' @param loss The loss function to use (\eqn{f(x)}). Supported values are: "mse".
#' @param prox_1 The proximal operator for the first function (\eqn{h(x)}).
#' @param prox_2 The proximal operator for the second function (\eqn{g(x)}).
#' @param pen_prox_1 The penalty for the first proximal operator. For the lasso, this would be \eqn{\lambda}. If operator does not include a penalty, set to 1.
#' @param pen_prox_2 The penalty for the second proximal operator. For the lasso, this would be \eqn{\lambda}. If operator does not include a penalty, set to 1.
#' @param prox_1_opts Optional arguments for first proximal operator. For the group lasso, this would be a vector of group IDs.
#' @param prox_2_opts Optional arguments for second proximal operator. For the group lasso, this would be a vector of group IDs.#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, as defined in Pedregosa et. al. (2018) as \eqn{\tau}.
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence threshold for the stopping criteria.
#' @param standardise Logical flag for standardising input data.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{x} \tab The solution to the OPT (See Pedregosa et. al. (2018)) \cr
#'    \tab \cr
#'    \code{u} \tab The solution to the dual problem (See Pedregosa et. al. (2018)) \cr
#'    \tab \cr
#'    \code{success} \tab Logical flag indicating whether ATOS converged, according to \code{tol}. \cr
#'    \tab \cr
#'    \code{num_it} \tab Number of iterations performed. If convergence is not reached, this will be \code{max_iter}. \cr
#'    \tab \cr
#' }
#'
#' @examples
#' # specify 6 groups of sizes 2, 3, and 4
#' group <- c(1, 1, 2, 2, 2, 3, 3, 3, 3,
#'            4, 4, 5, 5, 5, 6, 6, 6, 6)
#' # set the weight for each group to the square root of the group's size
#' wt <- rep(c(sqrt(2), sqrt(3), sqrt(4)), 2)
#' names(wt) <- 1:6
#' # compute different lambda sequences
#' lambda.max <- lambdaGroupSLOPE(method="max", fdr=0.1, group=group, wt=wt) 
#' lambda.mean <- lambdaGroupSLOPE(method="mean", fdr=0.1, group=group, wt=wt) 
#' lambda.corrected <- lambdaGroupSLOPE(method="corrected", fdr=0.1,
#'                                      group=group, wt=wt, n.obs=1000)
#' rbind(lambda.max, lambda.mean, lambda.corrected)
#' #                      [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#' # lambda.max       2.023449 1.844234 1.730818 1.645615 1.576359 1.517427
#' # lambda.mean      1.880540 1.723559 1.626517 1.554561 1.496603 1.447609
#' # lambda.corrected 1.880540 1.729811 1.637290 1.568971 1.514028 1.467551
#'
#' @references F. Pedregosa, G. Gidel (2018) \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

atos <- function(X, y, loss, prox_1, prox_2, pen_prox_1 = 0.5, pen_prox_2 = 0.5, max_iter = 5000, backtracking = 0.7, max_iter_backtracking = 100, tol = 1e-5,
                  prox_1_opts = NULL, prox_2_opts = NULL, group_ids, standardise=FALSE){

 if (standardise == TRUE){
   y = standardise_SGL(y) # standardise 
   X = center_scale(X, TRUE)$x # center-scale 
  }

  p = dim(X)[2]
  num_obs = dim(X)[1]
  x0 = rep(0,p) 
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  if (loss == "mse"){
    f = mse_loss
    f_grad = mse_grad
    f_opts=list(y=y,X=X,num_obs=num_obs)
    f_grad_opts=list(y=y,X=X,num_obs=num_obs)
  } else {print("Loss function not supported")}

  # initial step size
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)

  z = do.call(prox_1, c(list(x0, pen_prox_1*step_size), prox_1_opts))
  fz = do.call(f, c(list(z), f_opts))
  grad_fz = do.call(f_grad, c(list(z), f_grad_opts)) # gradient at z
  u = rep(0,p)
  x = do.call(prox_2, c(list(z - step_size * grad_fz, pen_prox_2*step_size), prox_2_opts))

  # Main fitting loop
  for (it in 1:max_iter){
    fz = do.call(f, c(list(z), f_opts))
    grad_fz = do.call(f_grad, c(list(z), f_grad_opts)) # gradient at z
    x = do.call(prox_2, c(list(z - step_size * (u + grad_fz), pen_prox_2*step_size), prox_2_opts))
    incr = x - z
    norm_incr = norm(incr,type="2")
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        x = do.call(prox_2, c(list(z - step_size * (u + grad_fz), pen_prox_2*step_size), prox_2_opts))
        incr = x - z
        norm_incr = norm(incr,type="2")
        rhs = fz +  crossprod(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol =  do.call(f, c(list(x), f_opts)) - rhs        
        if (ls_tol <= LS_EPS){
          break
        }
        else {
        step_size = step_size*backtracking 
        }
      }
    }

    z = do.call(prox_1, c(list(x + step_size * u, pen_prox_1*step_size), prox_1_opts)) 
    u = u + (x - z) / step_size
    certificate = norm_incr / step_size

    if (it > 0 & certificate < tol){ # Check for convergence
      success = 1
      break
    }
  }

  # Generate output 
  threshold_x = quantile(abs(x[which(x!=0)]))[4]*1e-5
  if (!is.na(threshold_x)) {
  threshold_x = ifelse(threshold_x<1e-5, 1e-5,0)
  x = ifelse(abs(x)>threshold_x,x,0)
  } 
  threshold_u = quantile(abs(u[which(u!=0)]))[4]*1e-5
  if (!is.na(threshold_u)) { 
  threshold_u = ifelse(threshold_u<1e-5, 1e-5,0)
  u = ifelse(abs(u)>threshold_u,u,0)
  } 
  out = c()
  out$x = as.matrix(x)
  out$num_iter = it
  out$success = success
  out$u = u
  class(out) <- "ATOS"
  return(out)
}
