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

#' print a `"sgs"` object
#'
#' Performs prediction from an [fit_sgs()] model fit.
#'
#' @param fit Object an object of class `"sgs"` from a call to [fit_sgs()] or [fit_sgs_cv()].
#' @param digits Number of digits to print in output.
#' @param ... ignored.
#' 
#' @seealso [fit_sgs()], [fit_sgs_cv()]
#' @family SGS-methods
#' 
#' @return A summary of the model fit. 
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
#' model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)
#' # print model
#' print(model)
#' @export

print.sgs <- function(fit, digits = max(3, getOption("digits") - 3), ...){ 
  num.nonzero <- apply(fit$beta,2, function(z){sum(z != 0)})
  cat("\n regression type: ", fit$type, "\n\n")
  print(cbind(lambdas = fit$lambdas, num.nonzero = num.nonzero, convergence = fit$success))
}

print.sgs_cv <- function(fit, digits = max(3, getOption("digits") - 3), ...){ 
  cat("\n regression type: ", fit$type, "\n\n")
  print(cbind(lambda = fit$errors$lambda, error = fit$errors$error_criteria, estimated_non_zero = fit$errors$num_non_zero))
}