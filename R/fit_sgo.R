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

#' Fit an SGO model.
#' 
#' Sparse-group OSCAR (SGO) main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#' 
#' \code{fit_sgo()} fits an SGO model (Feser and Evangelou (2024)) using adaptive three operator splitting (ATOS). SGO uses the same model set-up as for SGS, but with different weights (see Bao et al. (2020) and Feser and Evangelou (2024)).
#' The penalties are given by (for a group \eqn{g} and variable \eqn{i}, with \eqn{p} variables and \eqn{m} groups):
#' \deqn{
#'   v_i = \sigma_1 + \sigma_2(p-i), \; w_g = \sigma_1 + \sigma_3(m-g),
#' }
#' where
#' \deqn{
#'      \sigma_1 = d_i\|X^\intercal y\|_\infty, \; \sigma_2 = \sigma_1/p, \; \sigma_3 = \sigma_1/m, \; d_i = i \times \exp{(-2)}. 
#' }
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
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between OSCAR and gOSCAR. Must be between 0 and 1. Recommended value is 0.95.
#' @param max_iter Maximum number of ATOS iterations to perform. 
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence tolerance for the stopping criteria.
#' @param standardise Type of standardisation to perform on \code{X}: 
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/\sqrt{n}}.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/n}.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"noBaone"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param screen Logical flag for whether to apply screening rules (see Feser and Evangelou (2024)). Screening discards irrelevant groups before fitting, greatly improving speed.
#' @param verbose Logical flag for whether to print fitting information.
#' @param v_weights Optional vector for the variable penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' @param w_weights Optional vector for the group penalty weights. Overrides the OSCAR penalties when specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' 
#' @return A list containing:
#' \item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{z}, which is usually the former. A filter is applied to remove very small values, where ATOS has not been able to shrink exactly to zero. Check this against \code{x} and \code{z}.}
#' \item{group_effects}{The group values from the regression. Taken by applying the \eqn{\ell_2} norm within each group on \code{beta}.}
#' \item{selected_var}{A list containing the indicies of the active/selected variables for each \code{"lambda"} value. Index 1 corresponds to the first column in X.}
#' \item{selected_grp}{A list containing the indicies of the active/selected groups for each \code{"lambda"} value. Index 1 corresponds to the first group in the \code{groups} vector. You can see the group order by running \code{unique(groups)}.}
#' \item{num_it}{Number of iterations performed. If convergence is not reached, this will be \code{max_iter}.}
#' \item{success}{Logical flag indicating whether ATOS converged, according to \code{tol}.}
#' \item{certificate}{Final value of convergence criteria.}
#' \item{x}{The solution to the original problem (see Pedregosa and Gidel (2018)).}
#' \item{z}{The updated values from applying the first proximal operator (see Pedregosa and Gidel (2018)).}
#' \item{u}{The solution to the dual problem (see Pedregosa and Gidel (2018)).}
#' \item{screen_set_var}{List of variables that were kept after screening step for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{S}_v} in Feser and Evangelou (2024)).}
#' \item{screen_set_grp}{List of groups that were kept after screening step for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{S}_g} in Feser and Evangelou (2024)).}
#' \item{epsilon_set_var}{List of variables that were used for fitting after screening for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{E}_v} in Feser and Evangelou (2024)).}
#' \item{epsilon_set_grp}{List of groups that were used for fitting after screening for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{E}_g} in Feser and Evangelou (2024)).}
#' \item{kkt_violations_var}{List of variables that violated the KKT conditions each \code{"lambda"} value. (corresponds to \eqn{\mathcal{K}_v} in Feser and Evangelou (2024)).}
#' \item{kkt_violations_grp}{List of groups that violated the KKT conditions each \code{"lambda"} value. (corresponds to \eqn{\mathcal{K}_g} in Feser and Evangelou (2024)).}
#' \item{pen_slope}{Vector of the variable penalty sequence.}
#' \item{pen_gslope}{Vector of the group penalty sequence.}
#' \item{screen}{Logical flag indicating whether screening was performed.}
#' \item{type}{Indicates which type of regression was performed.}
#' \item{intercept}{Logical flag indicating whether an intercept was fit.}
#  \item{standardise}{Indicates the type of standardisation used.}
#' \item{lambda}{Value(s) of \eqn{\lambda} used to fit the model.}
#'
#' @family SGS-methods
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run SGO
#' model = fit_sgo(X = data$X, y = data$y, groups = groups, type="linear", path_length = 5, 
#' alpha=0.95, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' @references Bao, R., Gu B., Huang, H. (2020). \emph{Fast OSCAR and OWL Regression via Safe Screening Rules}, \url{https://proceedings.mlr.press/v119/bao20b}
#' @references Feser, F., Evangelou, M. (2023). \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
#' @references Pedregosa, F., Gidel, G. (2018). \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

fit_sgo <- function(X, y, groups, type="linear", lambda="path", path_length=20, min_frac=0.05, alpha=0.95, max_iter=5000, backtracking=0.7, max_iter_backtracking=100, tol=1e-5, standardise="l2", intercept=TRUE, screen=TRUE, verbose=FALSE, w_weights=NULL, v_weights=NULL){
  if (is.null(v_weights) & is.null(w_weights)){
    # create pen weights
    p = ncol(X)
    m = length(unique(groups))
    sigma_1 = exp(-2)*norm(t(X)%*%y,type="I")
    sigma_2 = sigma_1/p
    sigma_3 = sigma_1/m
    v_weights = sapply(1:p, function(x) sigma_1 + sigma_2*(p-x))
    w_weights = sapply(1:m, function(x) sigma_1 + sigma_3*(m-x))
  }

  # check group ordering to ensure no gaps in group numbers and ordered from 1 to m
  if (check_group_vector(groups)){
    reorder_id = FALSE
    ordered_grp_ids = groups
  } else {
    reorder_id = TRUE
    grp_new = reorder_group(groups)
    order_grp = order(grp_new,decreasing=FALSE)
    ordered_grp_ids = match(groups[order_grp], sort(unique(groups[order_grp])))
    X = X[,order_grp]
  }

  # Run main fitting function
  out = general_fit(X, y, ordered_grp_ids, "sgs", gen_path_sgs, sgs_var_screen, sgs_grp_screen, sgs_kkt_check, type, lambda, path_length, alpha, 0.1, 0.1, 3, 
                      backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, v_weights, w_weights, screen, verbose, FALSE, FALSE)
  
  # put group ordering back to original
  if (reorder_id){
    out = reorder_output(out, intercept, order_grp, groups)
  }

  return(out)
}
