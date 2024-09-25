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

#' Extracts coefficients for one of the following object types: `"sgs"`, `"sgs_cv"`, `"gslope"`, `"gslope_cv"`.
#'
#' Print the coefficients using model fitted with one of the following functions: [fit_sgs()], [fit_sgs_cv()], [fit_gslope()], [fit_gslope_cv()], [fit_sgo()], [fit_sgo_cv()], [fit_goscar()], [fit_goscar_cv()].
#' The predictions are calculated for each \code{"lambda"} value in the path.
#'
#' @param object Object of one of the following classes: \code{"sgs"}, \code{"sgs_cv"}, \code{"gslope"}, \code{"gslope_cv"}.
#' @param ... further arguments passed to stats function.
#' 
#' @seealso [fit_sgs()], [fit_sgs_cv()], [fit_gslope()], [fit_gslope_cv()]
#' @family SGS-methods
#' @family gSLOPE-methods
#' 
#' @return The fitted coefficients
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
#' model_coef = coef(model)

#' @method coef sgs
#' @export
coef.sgs <- function(object, ...){
    return(object$beta)
}

#' @method coef sgs_cv
#' @export
coef.sgs_cv <- function(object, ...){
    return(object$beta)
}

#' @method coef gslope
#' @export
coef.gslope <- function(object, ...){
    return(object$beta)
}

#' @method coef gslope_cv
#' @export
coef.gslope_cv <- function(object, ...){
    return(object$beta)
}