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

#' plot a `"sgs_cv"` object
#'
#' Plots the pathwise solution of a cross-validation fit, from a call to [fit_sgs_cv()] 
#'
#' @param x Object an object of class \code{"sgs_cv"} from a call to [fit_sgs()].
#' @param how_many Defines how many predictors to plot. Plots the predictors in decreasing order of largest absolute value.
#' @param ... further arguments passed to base function.
#' 
#' @seealso [fit_sgs_cv()]
#' @family SGS-methods
#' 
#' @return A list containing:
#' \item{response}{The predicted response. In the logistic case, this represents the predicted class probabilities.}
#' \item{class}{The predicted class assignments. Only returned if type = "logistic" in the \code{"sgs"} object.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,2,2,3)
#' # generate data
#' data = generate_toy_data(p=5, n=4, groups = groups, seed_id=3,signal_mean=20,group_sparsity=1)
#' # run SGS 
#' cv_model = fit_sgs_cv(X = data$X, y = data$y, groups=groups, type = "linear", 
#' nlambda = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, 
#' min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=FALSE)
#' plot(cv_model, how_many = 10)
#' @export

plot.sgs_cv <- function(x, how_many = 10, ...){ 
  beta_matrix = matrix(0,nrow=length(x$fit$z),ncol=length(x$errors$lambda))
  for (i in 1:length(x$errors$lambda)){
    if (x$fit$intercept){
      beta_matrix[,i] = as.vector(x$all_models[[i]]$beta[-1])
    } else {
      beta_matrix[,i] = as.vector(x$all_models[[i]]$beta)
    }
  }
  plot_path(beta_matrix=beta_matrix,lambdas=x$errors$lambda,how_many=how_many,main="Pathwise solution")
}