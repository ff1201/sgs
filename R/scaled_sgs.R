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

#' fits a scaled SGS model
#'
#' Fits an SGS model using the noise estimation procedure (Algorithm 5 from Bogdan). This estimates \eqn{\lambda} and then fits the model using the estimated value. It is an alternative approach to cross-validation ([fit_sgs_cv()]).
#'
#' @param X Input matrix of dimensions \eqn{p \times n}.
#' @param y Output vector of dimension \eqn{n}.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param type The type of regression to perform. Supported values are: "linear" and "logistic".
#' @param pen_method The type of penalty sequences to use.
#'   - `1` uses the vMean SGS and gMean gSLOPE sequences. 
#'   - `2` uses the vMax SGS and gMean gSLOPE sequences.
#'   - `3` uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param alpha The value of \eqn{\alpha}, defines the convex balance between SLOPE and gSLOPE.
#' @param vFDR Defines the desired variable FDR level, which determines the shape of the variable penalties.
#' @param gFDR Defines the desired group FDR level, which determines the shape of the group penalties.
#' @param standardise Type of standardisation. 
#'   - `l2` standardises the input data to have \eqn{\ell_2} norms of one.
#'   - `l1` standardises the input data to have \eqn{\ell_1} norms of one.
#'   - `sd` standardises the input data to have standard deviation of one.
#'   - `none` no standardisation applied.
#' @param intercept Logical flag for whether to fit an intercept.
#' @param verbose Logical flag for whether to print fitting information.
#' 
#' @return An object of type `"sgs"` containing model fit information.
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
#' # run noise estimation 
#' model = noise_estimation(X=data$X, y=data$y, groups=groups)
#' @references F. Feser, M. Evangelou \emph{Sparse-group SLOPE}, \url{https://github.com/ff1201/sgs}
#' @export

scaled_sgs <- function(X, y, groups, type="linear", pen_method = 1, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise="l2", intercept=TRUE, verbose=FALSE){
  # Run Algorithm 5 of Section 3.2.3. (Bogdan et al.) for SGS
      if (intercept) {
        selected <- 1
        X_2 = cbind(1,X)
        y_1 = y-mean(y)
      } else {
        selected <- integer(0)
      }

      out=standardise_sgs(X=X,y=y,standardise,intercept,dim(X)[1])
      attempts = 0
      repeat {
        selected_prev <- selected

        noise_est <- estimateNoise(X_2[, selected], y_1, intercept)
        
        fit <- fit_sgs(X=X, y=y, groups=groups, pen_method=pen_method, lambda=noise_est*out$scale_pen, alpha=alpha, vFDR=vFDR, gFDR=gFDR,intercept=intercept)
        
        attempts = attempts + 1
        if (verbose){print(paste0("Loop number: ", attempts))}

        selected <- fit$selected_var
        if (intercept) {
          selected <- union(1, selected)
        }

        if (identical(selected, selected_prev)) {
          break
        }

        if (length(selected) + 1 >= n) {
          stop("selected >= n-1 variables; cannot estimate variance")
        }
    }  
return(fit)
} 