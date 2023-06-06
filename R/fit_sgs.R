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

#' fit an SGS model
#' 
#' Sparse-group SLOPE (SGS) main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#' 
#' \code{fit_sgs()} fits an SGS model using adaptive three operator splitting (ATOS). SGS is a sparse-group method, so that it selects both variables and groups. Unlike group selection approaches, not every variable within a group is set as active.
#' It solves the convex optimisation problem given by 
#' \deqn{
#'   \frac{1}{2n} f(b ; y, \mathbf{X}) + \lambda \alpha \sum_{i=1}^{p}v_i |b|_{(i)} + \lambda (1-\alpha)\sum_{g=1}^{m}w_g \sqrt{p_g} \|b^{(g)}\|_2,
#' }
#' where \eqn{f(\cdot)} is the loss function. In the case of the linear model, the loss function is given by the mean-squared error loss:
#' \deqn{
#'  f(b; y, \mathbf{X}) = \left\|y-\mathbf{X}b \right\|_2^2.
#' }
#' In the logistic model, the loss function is given by 
#' \deqn{
#' f(b;y,\mathbf{X})=-1/n \log(\mathcal{L}(b; y, \mathbf{X})).
#' }
#' where the log-likelihood is given by
#' \deqn{
#'  \mathcal{L}(b; y, \mathbf{X}) = \sum_{i=1}^{n}\left\{y_i b^\intercal x_i - \log(1+\exp(b^\intercal x_i)) \right\}.
#' }
#' SGS can be seen to be a convex combination of SLOPE and gSLOPE, balanced through \code{alpha}, such that it reduces to SLOPE for \code{alpha = 0} and to gSLOPE for \code{alpha = 1}. 
#' The penalty parameters in SGS are sorted so that the largest coefficients are matched with the largest penalties, to reduce the FDR.
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param pen_method The type of penalty sequences to use (see Feser et al. (2023)):
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"3"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param lambda The value of \eqn{\lambda}, which defines the level of sparsity in the model. Can be picked using cross-validation (see [fit_sgs_cv()]). Must be a positive value.
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between SLOPE and gSLOPE. Must be between 0 and 1.
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties. Must be between 0 and 1.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
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
#' @param v_weights Optional vector for the variable penalty weights. Overrides the penalties from \code{pen_method} if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' @param w_weights Optional vector for the group penalty weights. Overrides the penalties from \code{pen_method} if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' @param x0 Optional initial vector for \eqn{x_0}.
#' @param u Optional initial vector for \eqn{u}.
#' @param verbose Logical flag for whether to print fitting information.
#'
#' @return A list containing:
#' \item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{u}, which is usually the former.}
#' \item{x}{The solution to the original problem (see Pedregosa et. al. (2018)).}
#' \item{u}{The solution to the dual problem (see Pedregosa et. al. (2018)).}
#' \item{z}{The updated values from applying the first proximal operator (see Pedregosa et. al. (2018)).}
#' \item{type}{Indicates which type of regression was performed.}
#' \item{pen_slope}{Vector of the variable penalty sequence.}
#' \item{pen_gslope}{Vector of the group penalty sequence.}
#' \item{lambda}{Value of \eqn{\lambda} used to fit the model.}
#' \item{success}{Logical flag indicating whether ATOS converged, according to \code{tol}.}
#' \item{num_it}{Number of iterations performed. If convergence is not reached, this will be \code{max_iter}.}
#' \item{certificate}{Final value of convergence criteria.}
#' \item{intercept}{Logical flag indicating whether an intercept was fit.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data = generate_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run SGS 
#' model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
#' vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' @references F. Feser, M. Evangelou \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references F. Pedregosa, G. Gidel (2018) \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

fit_sgs <- function(X, y, groups, pen_method = 1, type="linear", lambda, alpha=0.95, vFDR=0.1, gFDR=0.1, max_iter = 5000, backtracking = 0.7, max_iter_backtracking = 100, tol = 1e-5, standardise = "l2", intercept = TRUE,w_weights=NULL,v_weights=NULL,x0 = NULL, u = NULL,verbose=FALSE){
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
    stop("using a binary variable with linear regression. Use logistic instead")
  }
  if (lambda<0){
    stop("lambda can not be negative")
  }
  if (alpha<0 | alpha>1){
    stop("alpha must be in [0,1]")
  }
  if (vFDR<=0 | vFDR>=1 | gFDR<=0 | gFDR>=1){
    stop("FDR must be in (0,1)")
  }
  if (alpha == 1 & lambda != 0){
    warning("this package is not optimised for SLOPE. Consider using the SLOPE package instead")
  }
  if (alpha == 0 & lambda != 0){
    warning("this package is not optimised for group SLOPE")
  }
  if (lambda == 0){
    warning("this package is not optimised for OLS. Consider using the lm function instead")
  }

  # -------------------------------------------------------------
  # pre-process data
  # ------------------------------------------------------------- 
  num_vars = dim(X)[2]
  if (type == "logistic" & intercept){num_vars = dim(X)[2]+1} # fit intercept directly for logistic model
  num_obs = dim(X)[1]
  groups_pen = groups

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
  if (standardise=="none" & intercept==FALSE){
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
      if (type == "logistic" & intercept){ 
        X = cbind(1,X)
        groups = c(1,groups+1)
      }
  }

  # -------------------------------------------------------------
  # set values
  # ------------------------------------------------------------- 
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
  lambda_org = lambda
  lambda = scale_pen*lambda

  # weights
  if (is.null(v_weights) & is.null(w_weights)){
    pens_out = generate_penalties(gFDR, vFDR, pen_method, groups_pen, alpha)
    pen_slope_org = pens_out$pen_slope_org
    pen_gslope_org = pens_out$pen_gslope_org
    pen_slope = alpha*lambda*pen_slope_org
    pen_gslope = (1-alpha)*lambda*pen_gslope_org
  } else {
    pen_slope_org = v_weights
    pen_gslope_org = w_weights
    pen_slope = alpha*lambda*v_weights
    pen_gslope = (1-alpha)*lambda*w_weights
  } 

  # penalty checks
  if (any(pen_slope < 0) | any(pen_gslope < 0)){
    stop("penalty sequences must be positive")
  }

  if (!is.decreasing(pen_slope) | !is.decreasing(pen_slope)){
    stop("penalty sequences must be decreasing")
  }

  # -------------------------------------------------------------
  # fitting
  # ------------------------------------------------------------- 
  if (type == "logistic" & intercept){
    fit_out = run_atos_log_inter(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, x0, u,
                          pen_slope, pen_gslope, f, f_grad,f_opts, f_grad_opts, crossprod_mat,verbose,groups_pen)
  } else {  
    fit_out = run_atos(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, x0, u,
                          pen_slope, pen_gslope, f, f_grad,f_opts, f_grad_opts, crossprod_mat,verbose)
  }
  if (fit_out$success == 0){ # check for convergence
    warning("algorithm did not converge, try using more iterations")
  } else {if (verbose==TRUE){print("Algorithm converged")}}

  # -------------------------------------------------------------
  # generate output
  # ------------------------------------------------------------- 
  out = c()
  x = fit_out$x
  z = fit_out$z
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

  if (type == "logistic" & intercept){
    X = X[,-1]
    log_intercept = out$beta[1]
    out$beta = out$beta[-1]
    groups = groups[-1]-1
    num_vars = num_vars - 1
  }
  # scale beta depending on transformations
  if (standardise != "none"){ 
    out$beta = out$beta/X_scale
    out$x_beta = out$x_beta/X_scale
  }

  if (length(out$beta[out$beta!=0]) != 0){
    if (min(abs(out$beta[out$beta!=0])) < 1e-3){ # remove small resid values if solution not stable
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>(1e-4))]))[4]*1e-3
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>=threshold_x)]))[4]*1e-2
      if (!is.na(threshold_x) & lambda_org!=0) { # When lambda = 0, we don't want to remove small values, as no penalisation is occuring
        threshold_x = ifelse(threshold_x>1e-2,1e-2, threshold_x) # if threshold too big, set to 1e-2 
        out$beta = ifelse(abs(out$beta)>threshold_x,out$beta,0)
        out$beta = as.matrix(out$beta)
      }
    }
  }

  out$selected_var = which(out$beta!=0)
  which_groups_out = which_groups(beta=out$beta,groups=groups)
  out$selected_group = which_groups_out[[1]]
  out$group.effects = which_groups_out[[2]]

  if (intercept){ # get beta back to original scale
    if (type == "linear"){
    out$beta = as.matrix(c(y_mean - sum(X_center*out$beta),out$beta))
    out$x_beta = as.matrix(c(y_mean - sum(X_center*out$x_beta),out$x_beta))
    } else {
    out$beta = as.matrix(c(log_intercept,out$beta))
    out$x_beta = as.matrix(c(log_intercept,out$x_beta))
    }
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
  out$u = fit_out$u
  out$type = type
  out$pen_slope = pen_slope_org
  out$pen_gslope = pen_gslope_org
  out$lambda=lambda_org
  out$success = fit_out$success
  out$num_it = fit_out$it
  out$certificate = fit_out$certificate
  out$intercept = intercept
  class(out) <- "sgs"
  
  return(out)
}