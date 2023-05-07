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

#' generate toy data
#'
#' Generates different types of datasets, which can then be fitted using sparse-group SLOPE. 
#'
#' The data is generated under a Gaussian linear model. The generated data can be grouped and sparsity can be provided at both a group and/or variable level.
#'
#' @param p The number of input variables.
#' @param n The number of observations.
#' @param rho Correlation coefficient. Must be in range \eqn{[0,1]}.
#' @param seed_id Seed to be used to generate the data matrix \eqn{X}.
#' @param grouped A logical flag indicating whether grouped data is required.
#' @param groups If item{grouped=TRUE}, the grouping structure is required. Each input variable should have a group id.
#' @param noise_level Defines the level of noise (\eqn{sigma}) to be used in generating the response vector \eqn{y}.
#' @param group_sparsity Defines the level of group sparsity. Must be in the range \eqn{[0,1]}.
#' @param var_sparsity Defines the level of variable sparsity. Must be in the range \eqn{[0,1]}. If \code{grouped=TRUE}, this defines the level of sparsity within each group, not globally.
#' @param data_mean Defines the mean of input predictors.
#' @param data_sd Defines the standard deviation of the signal (\eqn{beta}).
#' @param orthogonal Logical flag as to whether the input matrix should be orthogonal.
#' @param signal_mean Defines the mean of the signal (\eqn{beta}).
#' @param signal_sd Defines the standard deviation of the signal (\eqn{beta}).
#' 
#' @return A list containing:
#' \item{y}{The response vector.}
#' \item{X}{The input matrix.}
#' \item{true_beta}{The true values of \eqn{beta} used to generate the response.}
#' \item{true_grp_id}{Indices of which groups are non-zero in item{true_beta}.}
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
#'
#' @export

generate_toy_data <- function(p, n, rho=0, seed_id=2, grouped=TRUE, groups, noise_level = 1, group_sparsity = 0.1, var_sparsity = 0.5,orthogonal = FALSE,
  data_mean = 0,data_sd = 1,signal_mean = 0,signal_sd = sqrt(10)){

  # Generates normally distributed toy datasets
  X=matrix(0,nrow=n,ncol=p)
  colnames(X) = 1:p
  for (j in 1:p){
    colnames(X)[j] = paste("v",j,sep="") # giving variable names to check the ordering at the end
  }
  out = c()

  if (grouped == TRUE){
    group2 <- paste0("grp", groups)
    group.id_org <- getGroupID(group2)
    group.length_org <- sapply(group.id_org, FUN=length)
    n_group <- length(group.id_org) # this extracts the total number of groups:
  
    for (group_id in 1:n_group){
      grp_id = group.id_org[[group_id]]
      set.seed((seed_id^2)*(13+group_id^2))
      X[,grp_id] = as.matrix(rnorm_multi(n=n,vars=length(grp_id),mu=data_mean,sd=data_sd,r=rho))
    }

    ind.relevant <- sort(sample(1:n_group, n_group*group_sparsity)) # 10% group sparsity -  indices of relevant groups

    true_beta <- rep(0, p) # pick true beta values (those with an effect)
    for (j in ind.relevant) {
      g_length_current = group.length_org[j]
      actives = which(rbern(n=g_length_current,prob=var_sparsity)==1)
      if (length(actives) == 0){ } else{
      true_beta[group.id_org[[j]]][actives] <- rnorm(length(actives),mean=signal_mean,sd=signal_sd)}
      }
    y <- X %*% true_beta + rnorm(n,mean=0,sd=noise_level) # calculate response
    true_ids = which(true_beta!=0)
    out$X = X
    out$true_beta = true_beta
    out$y = y
    out$true_var_id = true_ids
    out$true_grp_id = ind.relevant
  }

  else { # not grouped
    set.seed(seed_id)
    X = as.matrix(rnorm_multi(n=n,vars=p,mu=data_mean,sd=data_sd,r=rho))
    true_beta <- rep(0, p) # pick true beta values (those with an effect)
    true_ids = which(rbinom(n=p,size=1,prob=var_sparsity)==1)

    true_beta[true_ids] = rnorm(length(true_ids),mean=signal_mean,sd=signal_sd)

    y <- X %*% true_beta + rnorm(n,mean=0,sd=noise_level) # calculate response
    true_ids = which(true_beta!=0)
    out$true_beta = true_beta
    out$y = y
    out$X = X
    out$true_var_id = true_ids
  }
  return(out)  
}