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

#' Fit an SGS model using k-fold cross-validation.
#'
#' Function to fit a pathwise solution of sparse-group SLOPE (SGS) models using k-fold cross-validation. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#'
#' Fits SGS models under a pathwise solution using adaptive three operator splitting (ATOS), picking the 1se model as optimum. Warm starts are implemented.
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param lambda The regularisation parameter. Defines the level of sparsity in the model. A higher value leads to sparser models: 
#'   - \code{"path"} computes a path of regularisation parameters of length \code{"path_length"}. The path will begin just above the value at which the first predictor enters the model and will terminate at the value determined by \code{"min_frac"}.
#'   - User-specified single value or sequence. Internal scaling is applied based on the type of standardisation. The returned \code{"lambda"} value will be the original unscaled value(s).
#' @param path_length The number of \eqn{\lambda} values to fit the model for. If \code{"lambda"} is user-specified, this is ignored.
#' @param min_frac Defines the termination point of the pathwise solution, so that \eqn{\lambda_\text{min} = min_frac \cdot \lambda_\text{max}}.
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between SLOPE and gSLOPE. Must be between 0 and 1. Recommended value is 0.95.
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties. Must be between 0 and 1.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
#' @param pen_method The type of penalty sequences to use (see Feser and Evangelou (2023)):
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"3"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param nfolds The number of folds to use in cross-validation.
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
#' @param error_criteria The criteria used to discriminate between models along the path. Supported values are: \code{"mse"} (mean squared error) and \code{"mae"} (mean absolute error).
#' @param screen Logical flag for whether to apply screening rules (see Feser and Evangelou (2024)). Screening discards irrelevant groups before fitting, greatly improving speed.
#' @param verbose Logical flag for whether to print fitting information.
#' @param v_weights Optional vector for the variable penalty weights. Overrides the penalties from pen_method if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}
#' @param w_weights Optional vector for the group penalty weights. Overrides the penalties from pen_method if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}

#' @return A list containing:
#' \item{all_models}{A list of all the models fitted along the path.}
#' \item{fit}{The 1se chosen model, which is a \code{"sgs"} object type.}
#' \item{best_lambda}{The value of \eqn{\lambda} which generated the chosen model.}
#' \item{best_lambda_id}{The path index for the chosen model.}
#' \item{errors}{A table containing fitting information about the models on the path.}
#' \item{type}{Indicates which type of regression was performed.}
#' 
#' @seealso [fit_sgs()]
#' @family model-selection
#' @family SGS-methods
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run SGS with cross-validation (the proximal functions can be found in utils.R)
#' cv_model = fit_sgs_cv(X = data$X, y = data$y, groups=groups, type = "linear", 
#' path_length = 5, nfolds=5, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, 
#' standardise="l2",intercept=TRUE,verbose=TRUE)
#' @references Feser, F., Evangelou, M. (2023). \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
#' @export

fit_sgs_cv = function(X, y, groups, type = "linear", lambda="path", path_length = 20, min_frac = 0.05, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, pen_method=1, nfolds=10, backtracking = 0.7, max_iter = 5000, max_iter_backtracking = 100, tol = 1e-5, standardise= "l2", intercept = TRUE, error_criteria = "mse", screen=TRUE, verbose = FALSE, v_weights = NULL, w_weights = NULL){
  out = general_fit_cv(X, y, groups, "sgs", gen_path_sgs, type, lambda, path_length, nfolds, alpha, vFDR, gFDR, pen_method, 
                      backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, v_weights, w_weights, 
                      error_criteria, screen, verbose, FALSE, FALSE)
  return(out)
}