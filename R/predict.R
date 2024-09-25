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

#' Predict using one of the following object types: `"sgs"`, `"sgs_cv"`, `"gslope"`, `"gslope_cv"`.
#'
#' Performs prediction from one of the following fits: [fit_sgs()], [fit_sgs_cv()], [fit_gslope()], [fit_gslope_cv()], [fit_sgo()], [fit_sgo_cv()], [fit_goscar()], [fit_goscar_cv()]. The predictions are calculated for each \code{"lambda"} value in the path.
#'
#' @param object Object of one of the following classes: \code{"sgs"}, \code{"sgs_cv"}, \code{"gslope"}, \code{"gslope_cv"}.
#' @param x Input data to use for prediction.
#' @param ... further arguments passed to stats function.
#' 
#' @seealso [fit_sgs()], [fit_sgs_cv()], [fit_gslope()], [fit_gslope_cv()], [fit_sgo()], [fit_sgo_cv()], [fit_goscar()], [fit_goscar_cv()]
#' @family SGS-methods
#' @family gSLOPE-methods
#' 
#' @return A list containing:
#' \item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
#' \item{class}{The predicted class assignments. Only returned if type = "logistic" in the model object.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run SGS 
#' model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
#' vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' # use predict function
#' model_predictions = predict(model, x = data$X)

#' @method predict sgs
#' @export
predict.sgs <- function(object, x, ...){
  if (object$type=="linear"){
    if (object$intercept){
      predictions = apply(object$beta, 2, function(z) arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$beta, 2, function(z) arma_mv(x,as.vector(z)))
    }
  } else if (object$type == "logistic"){
    predictions = c()
    if (object$intercept){
      predictions$response = apply(object$beta, 2, function(z) sigmoid(arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$beta, 2, function(z) sigmoid(arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  } 
  return(predictions)
}

#' @method predict sgs_cv
#' @export
predict.sgs_cv <-  function(object, x, ...){
  if (object$fit$type=="linear"){
    if (object$fit$intercept){
      predictions = apply(object$all_models$beta, 2, function(z) arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$all_models$beta, 2, function(z) arma_mv(x,as.vector(z)))
    }
  } else if (object$fit$type == "logistic"){
    predictions = c()
    if (object$fit$intercept){
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  }
  return(predictions)
}

#' @method predict gslope
#' @export
predict.gslope <- function(object, x, ...){
  if (object$type=="linear"){
    if (object$intercept){
      predictions = apply(object$beta, 2, function(z) arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$beta, 2, function(z) arma_mv(x,as.vector(z)))
    }
  } else if (object$type == "logistic"){
    predictions = c()
    if (object$intercept){
      predictions$response = apply(object$beta, 2, function(z) sigmoid(arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$beta, 2, function(z) sigmoid(arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  } 
  return(predictions)
}

#' @method predict gslope_cv
#' @export
predict.gslope_cv <-  function(object, x, ...){
  if (object$fit$type=="linear"){
    if (object$fit$intercept){
      predictions = apply(object$all_models$beta, 2, function(z) arma_mv(Matrix::cbind2(1,x),as.vector(z)))
    } else {
      predictions = apply(object$all_models$beta, 2, function(z) arma_mv(x,as.vector(z)))
    }
  } else if (object$fit$type == "logistic"){
    predictions = c()
    if (object$fit$intercept){
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(arma_mv(Matrix::cbind2(1,x),as.vector(z))))
    } else {
      predictions$response = apply(object$all_models$beta, 2, function(z) sigmoid(arma_mv(x,as.vector(z))))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  }
  return(predictions)
}
