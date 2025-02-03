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

#' Fit an SGS model.
#' 
#' Sparse-group SLOPE (SGS) main fitting function. Supports both linear and logistic regression, both with dense and sparse matrix implementations.
#' 
#' \code{fit_sgs()} fits an SGS model (Feser and Evangelou (2023)) using adaptive three operator splitting (ATOS). SGS is a sparse-group method, so that it selects both variables and groups. Unlike group selection approaches, not every variable within a group is set as active.
#' It solves the convex optimisation problem given by 
#' \deqn{
#'   \frac{1}{2n} f(b ; y, \mathbf{X}) + \lambda \alpha \sum_{i=1}^{p}v_i |b|_{(i)} + \lambda (1-\alpha)\sum_{g=1}^{m}w_g \sqrt{p_g} \|b^{(g)}\|_2,
#' }
#' where \eqn{f(\cdot)} is the loss function and \eqn{p_g} are the group sizes. The penalty parameters in SGS are sorted so that the largest coefficients are matched with the largest penalties, to reduce the FDR. 
#' For the variables: \eqn{|\beta|_{(1)}\geq \ldots \geq |\beta|_{(p)}} and \eqn{v_1 \geq \ldots \geq v_p \geq 0}. 
#' For the groups: \eqn{\sqrt{p_1}\|\beta^{(1)}\|_2 \geq \ldots\geq \sqrt{p_m}\|\beta^{(m)}\|_2} and \eqn{w_1\geq \ldots \geq w_g \geq 0}.
#' In the case of the linear model, the loss function is given by the mean-squared error loss:
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
#' SGS can be seen to be a convex combination of SLOPE and gSLOPE, balanced through \code{alpha}, such that it reduces to SLOPE for \code{alpha = 1} and to gSLOPE for \code{alpha = 0}. 
#' The penalty parameters in SGS are sorted so that the largest coefficients are matched with the largest penalties, to reduce the FDR.
#' For the group penalties, see [fit_gslope()]. For the variable penalties, the vMean SGS sequence (\code{pen_method=1}) (Feser and Evangelou (2023)) is given by 
#' \deqn{
#' v_i^{\text{mean}} = \overline{F}_{\mathcal{N}}^{-1} \left( 1 - \frac{q_v i}{2p} \right), \; \text{where} \; \overline{F}_{\mathcal{N}}(x) := \frac{1}{m} \sum_{j=1}^{m} F_{\mathcal{N}} \left( \alpha x + \frac{1}{3} (1-\alpha) a_j w_j \right),\; i = 1,\ldots,p,
#' }
#' where \eqn{F_\mathcal{N}} is the cumulative distribution functions of a standard Gaussian distribution. The vMax SGS sequence (\code{pen_method=2}) (Feser and Evangelou (2023)) is given by
#' \deqn{
#' v_i^{\text{max}} = \max_{j=1,\dots,m} \left\{ \frac{1}{\alpha} F_{\mathcal{N}}^{-1} \left(1 - \frac{q_v i}{2p}\right) - \frac{1}{3\alpha}(1-\alpha) a_j w_j \right\},
#' }
#' The BH SLOPE sequence (\code{pen_method=3}) (Bogdan et al. (2015)) is given by 
#' \deqn{
#' v_i = z(1-i q_v/2p),
#' }
#' where \eqn{z} is the quantile function of a standard normal distribution.
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
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between SLOPE and gSLOPE. Must be between 0 and 1. Recommended value is 0.95.
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties. Must be between 0 and 1.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
#' @param pen_method The type of penalty sequences to use (see Feser and Evangelou (2023)):
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"3"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
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
#' @param v_weights Optional vector for the variable penalty weights. Overrides the penalties from \code{pen_method} if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' @param w_weights Optional vector for the group penalty weights. Overrides the penalties from \code{pen_method} if specified. When entering custom weights, these are multiplied internally by \eqn{\lambda} and \eqn{1-\alpha}. To void this behaviour, set \eqn{\lambda = 2} and \eqn{\alpha = 0.5}.
#' 
#' @return A list containing:
#' \item{beta}{The fitted values from the regression. Taken to be the more stable fit between \code{x} and \code{z}, which is usually the former. A filter is applied to remove very small values, where ATOS has not been able to shrink exactly to zero. Check this against \code{x} and \code{z}.}
#' \item{group_effects}{The group values from the regression. Taken by applying the \eqn{\ell_2} norm within each group on \code{beta}.}
#' \item{selected_var}{A list containing the indicies of the active/selected variables for each \code{"lambda"} value. Index 1 corresponds to the first column in X.}
#' \item{selected_grp}{A list containing the indicies of the active/selected groups for each \code{"lambda"} value. Index 1 corresponds to the first group in the \code{groups} vector.}
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
#' # run SGS 
#' model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", path_length = 5, 
#' alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' @references Bogdan, M., van den Berg, E., Sabatti, C., Candes, E. (2015). \emph{SLOPE - Adaptive variable selection via convex optimization}, \url{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-9/issue-3/SLOPEAdaptive-variable-selection-via-convex-optimization/10.1214/15-AOAS842.full}
#' @references Feser, F., Evangelou, M. (2023). \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references Feser, F., Evangelou, M. (2024). \emph{Strong screening rules for group-based SLOPE models}, \url{https://arxiv.org/abs/2405.15357}
#' @references Pedregosa, F., Gidel, G. (2018). \emph{Adaptive Three Operator Splitting}, \url{https://proceedings.mlr.press/v80/pedregosa18a.html}
#' @export

fit_sgs <- function(X, y, groups, type="linear", lambda="path", path_length=20, min_frac=0.05, alpha=0.95, vFDR=0.1, gFDR=0.1, pen_method=1, max_iter=5000, backtracking=0.7, max_iter_backtracking=100, tol=1e-5, standardise="l2", intercept=TRUE, screen=TRUE, verbose=FALSE, w_weights=NULL, v_weights=NULL){
  out = general_fit(X, y, groups, "sgs", gen_path_sgs, sgs_var_screen, sgs_grp_screen, sgs_kkt_check, type, lambda, path_length, alpha, vFDR, gFDR, pen_method, 
                      backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, v_weights, w_weights, screen, verbose, FALSE, FALSE)
  return(out)
}
