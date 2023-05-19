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

#' generate penalty sequences for SGS
#'
#' Generates variable and group penalties for SGS.
#'
#' The vMean and vMax SGS sequences are variable sequences derived specifically to give variable false discovery rate (FDR) control for SGS under orthogonal designs (see Feser et. al. (2023)).
#' The BH SLOPE sequence is derived in Bodgan et. al. (2015) and has links to the Benjamini-Hochberg critical values. The sequence provides variable FDR-control for SLOPE under orthogonal designs.
#' The gMean gSLOPE sequence is derived in Brzyski et. al. (2015) and provides group FDR-control for gSLOPE under orthogonal designs.
#'
#' @param vFDR Defines the desired variable false discovery rate (FDR) level, which determines the shape of the variable penalties.
#' @param gFDR Defines the desired group false discovery rate (FDR) level, which determines the shape of the group penalties.
#' @param pen_method The type of penalty sequences to use (see Feser et. al. (2023)):
#'   - \code{"1"} uses the vMean SGS and gMean gSLOPE sequences. 
#'   - \code{"2"} uses the vMax SGS and gMean gSLOPE sequences.
#'   - \code{"3"} uses the BH SLOPE and gMean gSLOPE sequences, also known as SGS Original.
#' @param groups A grouping structure for the input data. Should take the form of a vector of group indices.
#' @param alpha The value of \eqn{\alpha}, defines the convex balance between SLOPE and gSLOPE.
#' 
#' @return A list containing:
#' \item{pen_slope_org}{A vector of the variable penalty sequence.}
#' \item{pen_gslope_org}{A vector of the group penalty sequence.}
#'
#' @examples
#' # specify a grouping structure
#' groups = c(rep(1:20, each=3),
#'           rep(21:40, each=4),
#'           rep(41:60, each=5),
#'           rep(61:80, each=6),
#'           rep(81:100, each=7))
#' # generate sequences
#' sequences = generate_penalties(gFDR=0.1, vFDR=0.1, pen_method=1, groups=groups, alpha=0.5)
#' 
#' @references F. Feser, M. Evangelou \emph{Sparse-group SLOPE: adaptive bi-level selection with FDR-control}, \url{https://arxiv.org/abs/2305.09467}
#' @references M. Bogdan, E. Van den Berg, C. Sabatti, W. Su, E. Candes (2015) \emph{SLOPE — Adaptive variable selection via convex optimization}, \url{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-9/issue-3/SLOPEAdaptive-variable-selection-via-convex-optimization/10.1214/15-AOAS842.full}
#' @references D. Brzyski, W. Su, M. Bodgdan (2015) \emph{Group SLOPE - adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1511.09078}
#' @export

generate_penalties <- function(gFDR, vFDR, pen_method, groups, alpha){
  num_vars = length(groups)
  group_ids = getGroupID(groups) 
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  num_groups = length(unique(groups))

  if (pen_method == 1){ # SGS variable mean
    pen_gslope_org = lambdaChiOrtho(fdr=gFDR, n.group=num_groups, wt=wt_per_grp,
                           group.sizes=len_each_grp, method="mean")
    if (alpha==0){pen_slope_org = rep(0,num_vars)}else{
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=1,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="mean")}
  } else if (pen_method == 2){ # SGS variable max
    pen_gslope_org = lambdaChiOrtho(fdr=gFDR, n.group=num_groups, wt=wt_per_grp,
                           group.sizes=len_each_grp, method="mean")
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=1,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="max")
  } else if (pen_method == 3){ # SGS original
    pen_gslope_org = lambdaChiOrtho(fdr=gFDR, n.group=num_groups, wt=wt_per_grp,
                           group.sizes=len_each_grp, method="mean")
    pen_slope_org = BH_sequence(q=vFDR,p=num_vars)
  } else {stop("method choice not valid")}
out=c()
out$pen_slope_org = pen_slope_org
out$pen_gslope_org = pen_gslope_org
return(out)
}