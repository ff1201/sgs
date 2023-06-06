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

#' fits the adaptively scaled SGS model (AS-SGS)
#'
#' Fits an SGS model using the noise estimation procedure, termed adaptively scaled SGS (Algorithm 2 from Feser et al (2023)). 
#' This adaptively estimates \eqn{\lambda} and then fits the model using the estimated value. It is an alternative approach to 
#' cross-validation ([fit_sgs_cv()]). The approach is only compatible with the SGS penalties.
#'
#' @param X Input matrix of dimensions \eqn{n \times p}{n*p}. Can be a sparse matrix (using class \code{"sparseMatrix"} from the \code{Matrix} package).
#' @param y Output vector of dimension \eqn{n}. For \code{type="linear"} should be continuous and for \code{type="logistic"} should be a binary variable.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: \code{"linear"} and \code{"logistic"}.
#' @param pen_method The type of penalty sequences to use.
#'   - \code{"1"} uses the vMean and gMean SGS sequences.
#'   - \code{"2"} uses the vMax and gMax SGS sequences.
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
#' @references F. Feser, M. Evangelou \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @export

as_sgs <- function(X, y, groups, type="linear", pen_method = 2, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise="l2", intercept=TRUE, verbose=FALSE){
    num_obs=dim(X)[1]
    if (intercept) {
        selected <- 1
        X_2 = cbind(1,X)
        y_1 = y-mean(y)
    } else {
        selected <- integer(0)
    }
    out=standardise_sgs(X=X,y=y,standardise,intercept,dim(X)[1])
    selected_prev = 1000
    attempts = 0
    repeat {
        selected_prev_2 = selected_prev
        selected_prev <- selected
        
        
        noise_est <- estimateNoise(X_2[, selected], y_1, intercept)
        pens_out = generate_penalties_2(gFDR, vFDR, pen_method=pen_method,groups,alpha,lambda = noise_est)
        
        fit <- fit_sgs(X=X, y=y, groups=groups, pen_method=pen_method, type, lambda=noise_est*out$scale_pen, alpha=alpha, vFDR=vFDR, gFDR=gFDR,intercept=intercept,
                       v_weights=pens_out$pen_slope_org,w_weights=pens_out$pen_gslope_org,standardise=standardise)
        
        selected <- fit$selected_var
        if (intercept) {
            selected <- union(1, selected+1)
        }
        
        if (identical(selected, selected_prev) |identical(selected, selected_prev_2)) {
            break
        }
        if (length(selected) + 1 >= num_obs) {
            break
        }
        attempts = attempts+1
        if (verbose){print(paste0("Loop number: ", attempts))}

        if (attempts >= 100){
            break
        }
    }  
    return(fit)
} 