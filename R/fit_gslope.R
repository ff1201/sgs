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

#' Fit a gSLOPE model.
#'
#' Group SLOPE (gSLOPE) main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#'
#' \code{fit_gslope()} fits a gSLOPE model (Brzyski et al. (2019)) using adaptive three operator splitting (ATOS). gSLOPE is a sparse-group method, so that it selects both variables and groups. Unlike group selection approaches, not every variable within a group is set as active.
#' It solves the convex optimisation problem given by
#' \deqn{
#'   \frac{1}{2n} f(b ; y, \mathbf{X}) + \lambda \sum_{g=1}^{m}w_g \sqrt{p_g} \|b^{(g)}\|_2,
#' }
#' where the penalty sequences are sorted and \eqn{f(\cdot)} is the loss function. In the case of the linear model, the loss function is given by the mean-squared error loss:
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
#' The penalty parameters in gSLOPE are sorted so that the largest group effects are matched with the largest penalties, to reduce the group FDR.
#' The gMean sequence (\code{pen_method=1}) is given by
#' \deqn{
#'  w_i^\text{mean} = \overline{F}^{-1}_{\chi_{p_j}} (1-q_gi/m),  \; i = 1,\dots,m,
#' \text{where} \; \overline{F}_{\chi_{p_j}}(x):= \frac{1}{m}\sum_{j=1}^{m}F_{\chi_{p_j}}(\sqrt{p_j}x),
#' }
#' where \eqn{F_{\chi_{p_j}}} is the cumulative distribution function of a \eqn{\chi} distribution with \eqn{p_j} degrees of freedom. The gMax sequence (\code{pen_method=2}) is given by
#' \deqn{
#'  w_i^{\text{max}} = \max_{j=1,\ldots,m} \left\{ \frac{1}{\sqrt{p_j}} F^{-1}_{\chi_{p_j}} \left( 1 - \frac{q_g i}{m} \right) \right\},
#' }
#' where \eqn{F_{\chi_{p_j}}} is the cumulative distribution function of a \eqn{\chi} distribution with \eqn{p_j} degrees of freedom.
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
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
#' @param pen_method The type of penalty sequences to use (see Brzyski et al. (2019)):
#'   - \code{"1"} uses the gMean gSLOPE sequence.
#'   - \code{"2"} uses the gMax gSLOPE sequence.
#' @param max_iter Maximum number of ATOS iterations to perform.
#' @param backtracking The backtracking parameter, \eqn{\tau}, as defined in Pedregosa and Gidel (2018).
#' @param max_iter_backtracking Maximum number of backtracking line search iterations to perform per global iteration.
#' @param tol Convergence tolerance for the stopping criteria.
#' @param standardise Type of standardisation to perform on \code{X}:
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/\sqrt{n}}.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one. When using this \code{"lambda"} is scaled internally by \eqn{1/n}.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"none"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param screen Logical flag for whether to apply screening rules (see Feser and Evangelou (2024)). Screening discards irrelevant groups before fitting, greatly improving speed.
#' @param verbose Logical flag for whether to print fitting information.
#' @param w_weights Optional vector for the group penalty weights. Overrides the penalties from \code{pen_method} if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda}. To void this behaviour, set \eqn{\lambda = 1}.
#' @param warm_start Optional list for implementing warm starts. These values are used as initial values in the fitting algorithm. Need to supply \code{"x"} and \code{"u"} in the form \code{"list(warm_x, warm_u)"}. Not recommended for use with a path or CV fit as start from the null model by design.
#' @param warm_start Optional list for implementing warm starts. These values are used as initial values in the fitting algorithm. Need to supply \code{"x"} and \code{"u"} in the form \code{"list(warm_x, warm_u)"}. Not recommended for use with a path or CV fit as start from the null model by design.
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
#' \item{u}{The solution to the dual problem (see Pedregosa and Gidel (2018)).}
#' \item{z}{The updated values from applying the first proximal operator (see Pedregosa and Gidel (2018)).}
#' \item{screen_set}{List of groups that were kept after screening step for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{S}} in Feser and Evangelou (2024)).}
#' \item{epsilon_set}{List of groups that were used for fitting after screening for each \code{"lambda"} value. (corresponds to \eqn{\mathcal{E}} in Feser and Evangelou (2024)).}
#' \item{kkt_violations}{List of groups that violated the KKT conditions each \code{"lambda"} value. (corresponds to \eqn{\mathcal{K}} in Feser and Evangelou (2024)).}
#' \item{pen_gslope}{Vector of the group penalty sequence.}
#' \item{screen}{Logical flag indicating whether screening was applied.}
#' \item{type}{Indicates which type of regression was performed.}
#' \item{intercept}{Logical flag indicating whether an intercept was fit.}
#' \item{standardise}{Type of standardisation used.}
#' \item{lambda}{Value(s) of \eqn{\lambda} used to fit the model.}
#' 
#' @family gSLOPE-methods
#'
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,1,2,2,3,3,3,4,4)
#' # generate data
#' data =  gen_toy_data(p=10, n=5, groups = groups, seed_id=3,group_sparsity=1)
#' # run gSLOPE
#' model = fit_gslope(X = data$X, y = data$y, groups = groups, type="linear", path_length = 5,
#' gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' @references Brzyski, D., Gossmann, A., Su, W., Bodgan, M. (2019). \emph{Group SLOPE – Adaptive Selection of Groups of Predictors}, \url{https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1411269}
#' @references Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
#' @references Pedregosa, F., Gidel, G. (2018). \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

fit_gslope <- function(X, y, groups, type="linear", lambda="path", path_length=20, min_frac=0.05, gFDR=0.1, pen_method=1, max_iter=5000, backtracking=0.7, max_iter_backtracking=100, tol=1e-5, standardise="l2", intercept=TRUE, screen=TRUE, verbose=FALSE, w_weights=NULL, warm_start = NULL){
  if (pen_method == 1){
    pen_method_gslope = 3
  } else if (pen_method == 2){
    pen_method_gslope = 4
  } else {
    stop("pen_method choice not valid")
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
  out = general_fit(X, y, ordered_grp_ids, "gslope", gen_path_gslope, NULL, gslope_grp_screen, gslope_kkt_check, type, lambda, path_length, 0, 0.1, gFDR, pen_method_gslope,
                      backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, NULL, w_weights, screen, verbose, FALSE, FALSE, warm_start)
  
  # put group ordering back to original
  if (reorder_id){
    out = reorder_output(out, intercept, order_grp, groups)
  }
  
  return(out)
}
