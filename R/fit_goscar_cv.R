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

#' Fit a gOSCAR model using k-fold cross-validation.
#'
#' Function to fit a pathwise solution of group OSCAR (gOSCAR) models using k-fold cross-validation. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#'
#' Fits gOSCAR models under a pathwise solution using adaptive three operator splitting (ATOS), picking the 1se model as optimum. Warm starts are implemented.
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param lambda The regularisation parameter. Defines the level of sparsity in the model. A higher value leads to sparser models: 
#'   - \code{"path"} computes a path of regularisation parameters of length \code{"path_length"}. The path will begin just above the value at which the first predictor enters the model and will terminate at the value determined by \code{"min_frac"}.
#'   - User-specified single value or sequence. Internal scaling is applied based on the type of standardisation. The returned \code{"lambda"} value will be the original unscaled value(s).
#' @param path_length The number of \eqn{\lambda} values to fit the model for. If \code{"lambda"} is user-specified, this is ignored.
#' @param min_frac Smallest value of \eqn{\lambda} as a fraction of the maximum value. That is, the final \eqn{\lambda} will be \code{"min_frac"} of the first \eqn{\lambda} value.
#' @param nfolds The number of folds to use in cross-validation.
#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).
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
#' @param w_weights Optional vector for the group penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda}. To void this behaviour, set \eqn{\lambda = 1}.
#' @param warm_start Optional list for implementing warm starts. These values are used as initial values in the fitting algorithm. Need to supply \code{"x"} and \code{"u"} in the form \code{"list(warm_x, warm_u)"}. Not recommended for use with a path or CV fit as start from the null model by design.
#'
#' @return A list containing:
#' \item{errors}{A table containing fitting information about the models on the path.}
#' \item{all_models}{Fitting information for all models fit on the path, which is a \code{"gslope"} object type.}
#' \item{fit}{The 1se chosen model, which is a \code{"gslope"} object type.}
#' \item{best_lambda}{The value of \eqn{\lambda} which generated the chosen model.}
#' \item{best_lambda_id}{The path index for the chosen model.}
#'
#' @seealso [fit_goscar()]
#' @family gSLOPE-methods
#' @family model-selection
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run gOSCAR with cross-validation
#' cv_model = fit_goscar_cv(X = data$X, y = data$y, groups=groups, type = "linear", path_length = 5, 
#' nfolds=5, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=TRUE)
#' @references Bao, R., Gu B., Huang, H. (2020). \emph{Fast OSCAR and OWL Regression via Safe Screening Rules}, \url{https://proceedings.mlr.press/v119/bao20b}
#' @references Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
#' @export

fit_goscar_cv = function(X, y, groups, type = "linear", lambda = "path", path_length = 20, min_frac = 0.05, nfolds=10, backtracking = 0.7, max_iter = 5000, max_iter_backtracking = 100, tol = 1e-5, standardise= "l2", intercept = TRUE, error_criteria = "mse", screen=TRUE, verbose = FALSE, w_weights = NULL, warm_start = NULL){
  if (is.null(w_weights)){
    # create pen weights
    m = length(unique(groups))
    sigma_1 = exp(-2)*norm(t(X)%*%y,type="I")
    sigma_3 = sigma_1/m
    w_weights = sapply(1:m, function(x) sigma_1 + sigma_3*(m-x))
  }
  out = general_fit_cv(X, y, groups, "gslope", gen_path_gslope, type, lambda, path_length, nfolds, 0, 0.1, 0.1, 3, 
                      backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, NULL, w_weights, 
                      error_criteria, screen, verbose, FALSE, FALSE, warm_start)
  return(out)
}