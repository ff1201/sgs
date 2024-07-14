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

#' Fits a scaled SGS model.
#'
#' Fits an SGS model using the noise estimation procedure (Algorithm 5 from Bogdan et. al. (2015)). This estimates \eqn{\lambda} and then fits the model using the estimated value. It is an alternative approach to cross-validation ([fit_sgs_cv()]).
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param pen_method The type of penalty sequences to use.
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"1"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param alpha The value of \eqn{\alpha}, which defines the convex balance between SLOPE and gSLOPE. Must be between 0 and 1.
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties. Must be between 0 and 1.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties. Must be between 0 and 1.
#' @param standardise Type of standardisation to perform on \code{X}: 
#'   - \code{"l2"} standardises the input data to have \eqn{\ell_2} norms of one.
#'   - \code{"l1"} standardises the input data to have \eqn{\ell_1} norms of one.
#'   - \code{"sd"} standardises the input data to have standard deviation of one.
#'   - \code{"none"} no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param verbose Logical flag for whether to print fitting information.
#' 
#' @return An object of type \code{"sgs"} containing model fit information (see [fit_sgs()]). 
#'
#' @seealso [as_sgs()]
#' @family model-selection
#' @family SGS-methods
#' 
#' @examples
#' # specify a grouping structure
#' groups = c(1,1,2,2,3)
#' # generate data
#' data =  gen_toy_data(p=5, n=4, groups = groups, seed_id=3,
#' signal_mean=20,group_sparsity=1,var_sparsity=1)
#' # run noise estimation 
#' model = scaled_sgs(X=data$X, y=data$y, groups=groups, pen_method=1)
#' @references Bogdan, M., Van den Berg, E., Sabatti, C., Su, W., Candes, E. (2015). \emph{SLOPE — Adaptive variable selection via convex optimization}, \url{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-9/issue-3/SLOPEAdaptive-variable-selection-via-convex-optimization/10.1214/15-AOAS842.full}
#' @export

scaled_sgs <- function(X, y, groups, type="linear", pen_method = 1, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise="l2", intercept=TRUE, verbose=FALSE){
  if (intercept) {
    selected <- 1
    X_2 = Matrix::cbind2(1,X)
    y_1 = y-mean(y)
  } else {
    selected = integer(0)
  }
  num_obs=nrow(X)
  out=standardise_data(X=X,y=y,standardise,intercept,nrow(X))
  attempts = 0
  repeat {
    selected_prev = selected

    noise_est = estimateNoise(X_2[, selected], y_1, intercept)
    
    fit = fit_sgs(X=X, y=y, groups=groups, pen_method=pen_method, lambda=noise_est*out$scale_pen, alpha=alpha, vFDR=vFDR, gFDR=gFDR,intercept=intercept,screen=FALSE)
    
    attempts = attempts + 1
    if (verbose){print(paste0("Loop number: ", attempts))}

    if (intercept) {
      selected = union(1, selected)
    }

    if (identical(selected, selected_prev)) {
      break
    }

    if (length(selected) + 1 >= num_obs) {
      stop("selected >= n-1 variables; cannot estimate variance")
    }
  }  
  return(fit)
} 