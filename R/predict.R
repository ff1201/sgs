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

#' predict using a \code{"sgs"} object
#'
#' Performs prediction from an [fit_sgs()] model fit.
#'
#' @param fit Object an object of class \code{"sgs"} from a call to [fit_sgs()].
#' @param predict_X Input data to use for prediction.
#' @param ... further arguments passed to stats function.
#' 
#' @seealso [fit_sgs()]
#' @family SGS-methods
#' 
#' @return A list containing:
#' item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
#' item{class}{The predicted class assignments. Only returned if type = "logistic" in the \code{"sgs"} object.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(rep(1:20, each=3),
#'           rep(21:40, each=4),
#'           rep(41:60, each=5),
#'           rep(61:80, each=6),
#'           rep(81:100, each=7))
#' # generate data
#' data = generate_toy_data(p=500, n=400, groups = groups, seed_id=3)
#' # run SGS 
#' model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, 
#' vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' # use predict function
#' model_predictions = predict(model, predict_X = data$X)
#' @export
#' @method predict sgs

predict.sgs <- function(fit, predict_X, ...){
  if (fit$type=="linear"){
    if (fit$intercept){
    predictions = arma_mv(cbind(1,predict_X),as.vector(fit$beta))
    } else {
      predictions = arma_mv(predict_X,as.vector(fit$beta))
    }
  } else if (fit$type == "logistic"){
    predictions = c()
    if (fit$intercept){
      predictions$response = sigmoid(arma_mv(cbind(1,predict_X),as.vector(fit$beta)))
    } else {
      predictions$response = sigmoid(arma_mv(predict_X,as.vector(fit$beta)))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  } else {
       stop("not a valid type")
  }
  return(predictions)
}