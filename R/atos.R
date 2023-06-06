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

#' adaptive three operator splitting (ATOS)
#'
#' Function for fitting adaptive three operator splitting (ATOS) with general convex penalties. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#'
#' \code{atos()} solves convex minimization problems of the form
#' \deqn{
#'   f(x) + g(x) + h(x),
#' }
#' where \eqn{f} is convex and differentiable with \eqn{L_f}-Lipschitz gradient, and \eqn{g} and \eqn{h} are both convex.
#' The algorithm is not symmetrical, but usually the difference between variations are only small numerical values, which are filtered out.
#' However, both variations should be checked regardless, by looking at \code{x} and \code{u}. An example for the sparse-group lasso (SGL) is given. 
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package)
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} needs to be continuous and for \code{type="logistic"} needs to be a binary variable.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param prox_1 The proximal operator for the first function, \eqn{h(x)}.
#' @param prox_2 The proximal operator for the second function, \eqn{g(x)}.
#' @param pen_prox_1 The penalty for the first proximal operator. For the lasso, this would be the sparsity parameter, \eqn{\lambda}. If operator does not include a penalty, set to 1.
#' @param pen_prox_2 The penalty for the second proximal operator. 
#' @param prox_1_opts Optional argument for first proximal operator. For the group lasso, this would be the group IDs. Note: this must be inserted as a list.
#' @param prox_2_opts Optional argument for second proximal operator. 
#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa et. al. (2018).
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence tolerance for the stopping criteria.
#' @param standardise Type of standardisation to perform on \code{X}: 
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"none"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param x0 Optional initial vector for \eqn{x_0}.
#' @param u Optional initial vector for \eqn{u}.
#' @param verbose Logical flag for whether to print fitting information.
#'
#' @return An object of class \code{"atos"} containing:
#' \item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{u}, which is usually the former.}
#' \item{x}{The solution to the original problem (see Pedregosa et. al. (2018)).}
#' \item{u}{The solution to the dual problem (see Pedregosa et. al. (2018)).}
#' \item{z}{The updated values from applying the first proximal operator (see Pedregosa et. al. (2018)).}
#' \item{type}{Indicates which type of regression was performed.}
#' \item{success}{Logical flag indicating whether ATOS converged, according to \code{tol}.}
#' \item{num_it}{Number of iterations performed. If convergence is not reached, this will be \code{max_iter}.}
#' \item{certificate}{Final value of convergence criteria.}
#' \item{intercept}{Logical flag indicating whether an intercept was fit.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(rep(1:20, each=3),
#'           rep(21:40, each=4),
#'           rep(41:60, each=5),
#'           rep(61:80, each=6),
#'           rep(81:100, each=7))
#' # define proximal operators
#' L1_prox <- function(input, lambda){ # Lasso proximal operator
#'  out = sign(input) * pmax(0, abs(input) - lambda)
#'  return(out)
#' }
#' group_L1_prox = function(input,lambda,group_info){ 
#'  n_groups = length(unique(group_info))
#'  out = rep(0,length(input))
#'  for (i in 1:n_groups){
#'    grp_idx = which(group_info == unique(group_info)[i])
#'    if (lambda == 0 & norm(input[grp_idx],type="2") == 0){ # 0/0 = 0
#'      out[grp_idx] = 0
#'    } else {
#'      out[grp_idx] = max((1-(lambda/norm(input[grp_idx],type="2"))),0) * input[grp_idx]}
#'  }
#'  return(out)
#' }
#' # generate data
#' data = generate_toy_data(p=500, n=400, groups = groups, seed_id=3)
#' # run atos (the proximal functions can be found in utils.R)
#' out = atos(X=data$X, y=data$y, type="linear", prox_1 = L1_prox, prox_2 = group_L1_prox, 
#' standardise="none", intercept=FALSE, prox_2_opts = list(groups))
#' @references F. Pedregosa, G. Gidel (2018) \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

atos <- function(X, y, type = "linear", prox_1, prox_2, pen_prox_1 = 0.5, pen_prox_2 = 0.5, max_iter = 5000, backtracking = 0.7, max_iter_backtracking = 100, tol = 1e-5,
                  prox_1_opts = NULL, prox_2_opts = NULL, standardise = "l2", intercept = TRUE,x0 = NULL, u = NULL,verbose=FALSE){
  # -------------------------------------------------------------
  # checks
  # ------------------------------------------------------------- 
  if (anyNA(y) | anyNA(X)) {
    stop("input contains missing values")
  }
  if (!(is.matrix(X) | is(X, 'sparseMatrix')) | !(is.matrix(y) | is.vector(y))){
    stop("X and y must be matrices/vectors. Use the as.matrix function to convert")
  }
  if (length(y) == 0) {
    stop("y is empty")
  }
  if (nrow(X) == 0) {
    stop("X is empty")
  }
  if (length(y) != nrow(X)) {
    stop("the number of samples in y must match the number of rows in X")
  }
  if (type == "logistic" & !is.binary(y)){
    stop("when using logistic regression the response, y, needs to be a binary variable")
  }
  if (type == "linear" & is.binary(y)){
    stop("using a binary variable with linear regression. use logistic instead")
  }
  # -------------------------------------------------------------
  # pre-process data
  # ------------------------------------------------------------- 
  num_vars = dim(X)[2]
  num_obs = dim(X)[1]

  if (sum(X==0) > (num_vars*num_obs)/2){
    warnings("X appears to be a sparse matrix. Try converting to dgCMatrix type for improved performance")
  }
  
  is_sparse = FALSE
  crossprod_mat = base::crossprod
  if (inherits(X,"dgCMatrix")){ # check if matrix is sparse
    crossprod_mat = Matrix::crossprod
    is_sparse = TRUE
  }

  if (standardise!="none" & is_sparse){
    stop("standardising a matrix that is sparse. this would remove sparsity of X. set standardise to none")
  }

  # standardise
  if (standardise=="none"){
      scale_pen = 1
      y_mean = 0 
      X_center = 0
      X_scale = 1
    } else {
      standardise_out = standardise_sgs(X=X,y=y,standardise=standardise,intercept=intercept,num_obs=num_obs,type=type)
      X = standardise_out$X
      X_scale = standardise_out$X_scale
      X_center = standardise_out$X_center
      y = standardise_out$y
      y_mean = standardise_out$y_mean
      scale_pen = standardise_out$scale_pen
  }

  # -------------------------------------------------------------
  # set values
  # ------------------------------------------------------------- 
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  # type of model
  if (type == "linear"){ 
    if (is_sparse){
      f = mse_loss_sparse
      f_grad = mse_grad_sparse
    } else {
      f = mse_loss
      f_grad = mse_grad
    }
  } else if (type == "logistic"){
    if (is_sparse){
      f = log_loss_sparse
      f_grad = log_grad_sparse
    } else {
      f = log_loss
      f_grad = log_grad
    }
  } else {stop("loss function not supported")} 
  f_opts = list(y=y,X=X,num_obs=num_obs)
  f_grad_opts = list(y=y,X=X,num_obs=num_obs)
  pen_prox_1 = scale_pen*pen_prox_1
  pen_prox_2 = scale_pen*pen_prox_2

  # -------------------------------------------------------------
  # initial fitting values
  # ------------------------------------------------------------- 
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
  z = do.call(prox_1, c(list(x0, pen_prox_1*step_size), prox_1_opts))
  fz = do.call(f, c(list(z), f_opts))
  grad_fz = do.call(f_grad, c(list(z), f_grad_opts)) # gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  x = do.call(prox_2, c(list(z - step_size * grad_fz, pen_prox_2*step_size), prox_2_opts))

  # -------------------------------------------------------------
  # fitting
  # ------------------------------------------------------------- 
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

    if (certificate < tol){ # Check for convergence
      success = 1
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }

  if (success == 0){ # check for convergence
    warning("algorithm did not converge, try using more iterations")
  } else {if (verbose==TRUE){print("Algorithm converged")}}

  # -------------------------------------------------------------
  # generate output
  # ------------------------------------------------------------- 
  out = c()
  out$x_beta = x
  if (max((x-z)^2) < 1e-3 & mean((x-z)^2) < 1e-3){ # if solutions are very similar, pick more stable version
    if (length(which(x!=0)) <= length(which(z!=0))){ # Picking the solution with less residual values, if this is true, x is picked
      out$beta = as.matrix(x)
    } else {
      out$beta = as.matrix(z)
    }
  } else { # if solutions aren't similar, pick x
    out$beta = as.matrix(x)
  }

  # scale beta depending on transformations
  if (standardise!="none"){ 
    out$beta = out$beta/X_scale
    out$x_beta = out$x_beta/X_scale
  }

  if (length(out$beta[out$beta!=0]) != 0){
    if (min(abs(out$beta[out$beta!=0])) < 1e-3){ # remove small resid values if solution not stable
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>(1e-4))]))[4]*1e-3
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>=threshold_x)]))[4]*1e-2
      if (!is.na(threshold_x) & !(pen_prox_1==0 & pen_prox_2==0)){ # When lambda = 0, we don't want to remove small values, as no penalisation is occuring
        threshold_x = ifelse(threshold_x>1e-2,1e-2, threshold_x) # if threshold too big, set to 1e-2 
        out$beta = ifelse(abs(out$beta)>threshold_x,out$beta,0)
        out$beta = as.matrix(out$beta)
      }
    }
  }

  out$selected_var = which(out$beta!=0)

  if (intercept){ # get beta back to original scale
    out$beta = as.matrix(c(y_mean - sum(X_center*out$beta),out$beta))
    out$x_beta = as.matrix(c(y_mean - sum(X_center*out$x_beta),out$x_beta))
  } 

  if (is.null(colnames(X))){ # Add variable names to output
    if (intercept){
        rownames(out$beta) = c("(Intercept)", paste0("v", 1:(num_vars)))
      } else {
        rownames(out$beta) = paste0("v", 1:num_vars)
      }
    } else {
    if (intercept){
        rownames(out$beta) = c("(Intercept)", colnames(X))
      } else {
        rownames(out$beta) = colnames(X)
      }
  }
  out$beta = as(out$beta,"CsparseMatrix")
  out$z = z
  out$x = x
  out$u = u
  out$type = type
  out$success = success
  out$num_it = it
  out$certificate = certificate
  out$intercept = intercept
  class(out) <- "atos"
  return(out)
}