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

#' fit an SGS model using CV
#'
#' Function to fit a pathwise solution of sparse-group SLOPE (SGS) models using k-fold cross-validation. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#'
#' Fits SGS models under a pathwise solution using adaptive three operator splitting (ATOS), picking the 1se model as optimum. Warm starts are implemented.
#'
#' @param X Input matrix of dimensions \eqn{p x n}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param pen_method The type of penalty sequences to use (see Feser et. al. (2023)):
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"3"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param nlambda The number of pathwise \eqn{\lambda} values to fit.
#' @param nfolds The number of folds to use in cross-validation.
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between SLOPE and gSLOPE. Must be between 0 and 1.
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties. Must be between 0 and 1.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa et. al. (2018).
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence tolerance for the stopping criteria.
#' @param min_frac Defines the termination point of the pathwise solution, so that \eqn{\lambda_\text{min} = min_frac \cdot \lambda_\text{max}}.
#' @param standardise Type of standardisation to perform on \code{X}: 
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"none"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param verbose Logical flag for whether to print fitting information.
#' @param v_weights Optional vector for the variable penalty weights. Overrides the penalties from pen_method if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}
#' @param w_weights Optional vector for the group penalty weights. Overrides the penalties from pen_method if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}
#' @param error_criteria The criteria used to discriminate between models along the path. Supported values are: \code{"mse"} (mean squared error) and \code{"mae"} (mean absolute error).
#' @param max_lambda Optional parameter, \eqn{\lambda_\text{max}}, which is used to fit the first model on the path. If not specificed, it is chosen to be just above the value which lets in the first variable (so that it is the null model). 
#'
#' @return A list containing:
#' \item{all_models}{A list of all the models fitted along the path.}
#' \item{fit}{The 1se chosen model, which is a \code{"sgs"} object type.}
#' \item{best_lambda}{The value of \eqn{\lambda} which generated the chosen model.}
#' \item{best_lambda_id}{The path index for the chosen model.}
#' \item{errors}{A table containing fitting information about the models on the path.}
#' \item{type}{Indicates which type of regression was performed.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data = generate_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run SGS with cross-validation (the proximal functions can be found in utils.R)
#' cv_model = fit_sgs_cv(X = data$X, y = data$y, groups=groups, type = "linear", 
#' nlambda = 5, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, 
#' standardise="l2",intercept=TRUE,verbose=TRUE)
#' @references F. Feser, M. Evangelou \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references F. Pedregosa, G. Gidel (2018) \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

fit_sgs_cv = function(X,y,groups,pen_method=1, type = "linear", nlambda = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, backtracking = 0.7, max_iter = 5000, max_iter_backtracking = 100, tol = 1e-5, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=FALSE, v_weights=NULL,w_weights=NULL,error_criteria = "mse",max_lambda = NULL){
  num_vars = ncol(X)
  num_obs = nrow(X)
  num_groups = length(unique(groups))
  group_ids = getGroupID(groups) 
  len_each_grp = sapply(group_ids, length)
  
  # weights
  if (is.null(v_weights) & is.null(w_weights)){
    pens_out = generate_penalties(gFDR, vFDR, pen_method, groups, alpha)
    pen_slope_org = pens_out$pen_slope_org
    pen_gslope_org = pens_out$pen_gslope_org
  } else {
    pen_slope_org = v_weights
    pen_gslope_org = w_weights
  }

  # standardise data
  standardise_out = standardise_sgs(X=X,y=y,standardise=standardise,intercept=intercept,num_obs=num_obs)
  X_path = standardise_out$X
  y_path = standardise_out$y
  
  # calculate path
  if (is.null(max_lambda)){
    lambdas = (1/standardise_out$scale_pen)*generate_lambda_path(X=as.matrix(X_path),y=as.matrix(y_path),groups=groups,alpha=alpha,min_frac=min_frac,nlambda=nlambda, v_weights=pen_slope_org, w_weights=pen_gslope_org,group.sizes = len_each_grp)
  } else {
    min_lambda = min_frac*max_lambda
    lambdas = exp(seq(log(max_lambda),log(min_lambda), (log(min_lambda) - log(max_lambda))/(nlambda-1))) 
  }

  # initialise CV variable for storing results
  all_data = data.frame(y,X)
  set.seed(5)
  folds = createFolds(y, k = nfolds, list=TRUE)
  all_errors = matrix(0,nrow=nlambda,ncol=nfolds)
  output_errors = data.frame(lambda=lambdas,error_criteria=rep(0,nlambda), num_non_zero = rep(0,nlambda))
  ## warm starts
  if (type=="logistic" & intercept){    
    initial_x = rep(0,num_vars+1)
    initial_u = rep(0,num_vars+1)
  } else {  
    initial_x = rep(0,num_vars)
    initial_u = rep(0,num_vars)}

  list_of_models = list()
  
  for (lambda_id in 1:nlambda){
    num_its = 0
    for (fold_id in 1:nfolds){
    
      # Create training design matrix and target data, leaving one out each time
      Train = all_data[as.numeric(unlist(folds[-fold_id])),]
      Train_y = Train$y
      Train_X = Train[,-1]

      # Create testing design matrix and target data
      Test = all_data[as.numeric(unlist(folds[fold_id])),]
      Test_y = Test$y
      Test_X = X[as.numeric(unlist(folds[fold_id])),]
    
      # Fit Model
      model = fit_sgs(X=as.matrix(Train_X), y=as.matrix(Train_y),groups=groups, type=type, lambda=lambdas[lambda_id], alpha=alpha, max_iter = max_iter, backtracking = backtracking, max_iter_backtracking = max_iter_backtracking, tol = tol, standardise=standardise, intercept=intercept,
        x0 = initial_x, u = initial_u, v_weights= pen_slope_org, w_weights = pen_gslope_org)

      # Error
      if (type=="linear"){
      if (intercept){
        if (error_criteria == "mse"){
          error_val = sum((Test_y-arma_mv(cbind(1,Test_X),as.vector(model$beta)))^2)}
          else if (error_criteria == "mae") {error_val = sum(abs(Test_y-arma_mv(cbind(1,Test_X),as.vector(model$beta))))} else{stop("not a valid criteria")}
      } else {
        if (error_criteria == "mse"){
          error_val = sum((Test_y-arma_mv(Test_X,as.vector(model$beta)))^2)
          } else if (error_criteria =="mae"){error_val = sum(abs(Test_y-arma_mv(Test_X,as.vector(model$beta))))} else {stop("not a valid criteria")}
      }
      } else if (type=="logistic"){
        if (intercept){
          error_val = 1-sum(ifelse(sigmoid(arma_mv(cbind(1,Test_X),as.vector(model$beta)))>=0.5,1,0) == Test_y)/length(Test_y)
        } else {
          error_val = 1-sum(ifelse(sigmoid(arma_mv(Test_X,as.vector(model$beta)))>=0.5,1,0) == Test_y)/length(Test_y)
        }
      }
      all_errors[lambda_id, fold_id] = error_val
      num_its = num_its + model$num_it
  }
      lambda_model = fit_sgs(X=X, y=y,groups=groups, type=type, lambda=lambdas[lambda_id], alpha=alpha, max_iter = max_iter, backtracking = backtracking, max_iter_backtracking = max_iter_backtracking, tol = tol, standardise=standardise, intercept=intercept,
      x0 = initial_x, u = initial_u, v_weights= pen_slope_org, w_weights = pen_gslope_org)
      list_of_models[[lambda_id]] = lambda_model
      # warm starts
      initial_x = lambda_model$x
      initial_u = lambda_model$u
      output_errors$error_criteria[lambda_id] = mean(all_errors[lambda_id,])
      output_errors$num_non_zero[lambda_id] = length(lambda_model$selected_var)
      if (verbose == TRUE){
        if (type == "linear"){print(paste0("Lambda ", lambda_id,"/",nlambda, " done. Lambda: ", round(lambdas[lambda_id],4), ". Number of non-zero: ",length(lambda_model$selected_var),". Error: ", output_errors$error_criteria[lambda_id], ". Avg iter: ", floor(num_its/nfolds)))}
        else if (type == "logistic"){
          print(paste0("Lambda ", lambda_id,"/",nlambda, " done. Lambda: ", round(lambdas[lambda_id],4), ". Number of non-zero: ",length(lambda_model$selected_var),". Misclassification error: ", output_errors$error_criteria[lambda_id], ". Avg iter: ", floor(num_its/nfolds)))
          } }
  }

  # Pick best lambda - 1se
  error_se = apply(all_errors,1,sd)/sqrt(nfolds)
  error_se_lambda_min = error_se[which.min(output_errors$error_criteria)]
  best_lambda = max(lambdas[output_errors$error_criteria < min(output_errors$error_criteria) + error_se_lambda_min])
 
  model = list_of_models[[match(best_lambda, lambdas)]]
  out = c()
  out$all_models = list_of_models
  out$fit = model
  out$best_lambda = best_lambda
  out$best_lambda_id = match(best_lambda, lambdas)
  out$errors = output_errors
  out$type = type
  class(out) <- "sgs_cv"
  return(out)
}