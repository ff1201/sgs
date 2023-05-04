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

###################################################
#                 Packages
###################################################
library(Rcpp)
library(Matrix)
library(SLOPE)
library(grpSLOPE)
library(MASS)
library(VGAM)
#library(RColorBrewer)
library(caret)
library(Rlab)
library(faux)
library(glmnet)
#library(pracma)

###################################################
#                 General functions
###################################################
## Stop 1 dim matrices becoming vectors
  old <- `[`
`[` <- function(...) { old(..., drop=FALSE) }
 
# matrix multiplcation
arma_code <-
  "arma::vec arma_mm(const arma::mat& m, const arma::vec& v) {
       return m * v;
   };"
arma_mm = cppFunction(code = arma_code, depends = "RcppArmadillo")

is.decreasing <- function(x){ # checks if a sequence is decreasing
  state = TRUE
  for (i in 2:length(x)){
    if (x[i] > x[i-1]){
      state = FALSE
    }
  }
  return(state)
} 

sigmoid = function(x) {
   1 / (1 + exp(-x))
}

BH_sequence = function(q,p){
  # Calculates the BH sequences for SLOPE and gSLOPE
  p_sequence = rep(0,p)
  for (i in 1:p){
    q_i = (i*q)/(2*p)
    p_sequence[i] = qnorm(1-q_i)
  }
return(p_sequence)
}

mse_grad <- function(y, X, input,num_obs){ # gradient of loss
  r = arma_mm(X,input) - y
  out = crossprod(X, r)/num_obs
  return(out)
}

mse_loss <- function(y, X, input, num_obs){ # loss function
  out = as.double(crossprod(y-arma_mm(X,input))/(2*num_obs))
  return(out)
}

mse_grad_sparse <- function(y, X, input,num_obs){ # gradient of loss
  r = X%*%input - y
  out = Matrix::crossprod(X, r)/num_obs
  return(out)
}

mse_loss_sparse <- function(y, X, input, num_obs){ # loss function
  out = as.double(Matrix::crossprod(y-X%*%input)/(2*num_obs))
  return(out)
}

#stable log implementations - from https://fa.bianp.net/blog/2019/evaluate_logistic/
logsig <- function(input){
  out = rep(0,length(input))
  idx_1 = input < -33
  idx_2 = input >= -33 & input < -18
  idx_3 = input >= -18 & input < 37
  idx_4 = input >= 37

  out[idx_1] = input[idx_1]
  out[idx_2] = input[idx_2] - exp(input[idx_2])
  out[idx_3] = -log1p(exp(-input[idx_3]))
  out[idx_4] = -exp(-input[idx_4])

  return(out)
}

log_loss <- function(y,X,input,num_obs){# stable version for y{0,1}
  #eps = 1e-15
  z = arma_mm(X,input)
  #z <- pmax(pmin(z, 1 - eps), eps)
  out = mean((1-y)*z - logsig(z))
  return(out)
}

expit_b <- function(t,b){
  out = rep(0,length(t))
  idx = t<0
  b_pos = b[idx]
  b_neg = b[!idx]
  exp_pos = exp(t[idx])
  exp_neg = exp(-t[!idx])
  out[idx] = ((1-b_pos)*exp_pos - b_pos)/(1+exp_pos)
  out[!idx] = ((1-b_neg) - b_neg*exp_neg)/(1+exp_neg)
  return(out)
}

log_grad <- function(y,X,input,num_obs){# stable version for y{0,1}
  z = arma_mm(X,input)
  s = expit_b(z,y)
  out = (arma_mm(t(X),s))/dim(X)[1]
  return(out)
}

log_loss_sparse <- function(y,X,input,num_obs) {
  out = log(1+exp(Matrix::crossprod(t(y),Matrix::crossprod(X,input))))/num_obs
  return(out)
}

log_grad_sparse <- function(y, X, input,num_obs) {
  out = (1/num_obs)*(-y/(1+exp(Matrix::crossprod(t(y),Matrix::crossprod(X,input)))))
  return(out)
}

init_lipschitz <- function(f, f_grad, x0, f_opts, f_grad_opts){
  L0 = 1e-3
  f0 = do.call(f, c(list(x0), f_opts))
  grad0 = do.call(f_grad, c(list(x0), f_grad_opts))

  x_tilde = x0 - (1 / L0)*grad0
  f_tilde = do.call(f, c(list(x_tilde), f_opts))

  for (i in 1:100){
    if (f_tilde <= f0){
      break
    } else {
      L0 = L0 * 10
      x_tilde = x0 - (1 / L0) * grad0
      f_tilde = do.call(f, c(list(x_tilde), f_opts))
    }
  }
  return(L0)
}

norm_vec <- function(x) sqrt(sum(x^2))

proxGroupSortedL1 <- function(y, lambda,group, group_id, ...) { 
  # proximal operator for group SLOPE - adapted so that the 0/0 = NaN error doesn't occur
  n.group = length(unique(group))

  if (length(lambda) != n.group) {
    stop("Length of lambda should be equal to the number of groups.")
  }

  # compute Euclidean norms for groups in y
  group.norm <- rep(NA, n.group)
  for (i in 1:n.group){
    selected <- group_id[[i]]
    group.norm[i] <- norm_vec(y[selected])
  }

  # get Euclidean norms of the solution vector
  prox.norm <- prox_sorted_L1(x=group.norm, lambda=lambda, ...)

  # compute the solution
  prox.solution <- rep(NA, length(y))
  for (i in 1:n.group){
    selected <- group_id[[i]]
    if (group.norm[i] == 0){ # to stop 0/0 = NaN
      prox.solution[selected] <- 0
      } else {
    prox.solution[selected] <- prox.norm[i] / group.norm[i] * y[selected] }
  }
  return(prox.solution)
}

sgs_convex_opt = function(X,y,beta,groups,num_obs,gslope_seq,slope_seq,intercept=TRUE){  
  # function to evaluation convex optimisation function
  beta_org = beta
  if (intercept==TRUE){beta = beta[-1]}
  num_groups = length(unique(groups))
  group_norm_sgs = rep(0,num_groups)
  group_lengths = rep(0,num_groups)
  len_so_far = 0
  for (i in 1:length(unique(groups))){
    g_length = table(groups)[i]
    group_lengths[i] = g_length
    group_norm_sgs[i] = norm(beta[(len_so_far+1):(len_so_far+g_length)],type="2")
    len_so_far = len_so_far+g_length
  }
  if (intercept==TRUE){
    loss = mse_loss(y=y,X=cbind(1,X),input=beta_org,num_obs = num_obs)
  } else {loss = mse_loss(y=y,X=X,input=beta_org,num_obs = num_obs)}
  group_norm_sgs = sqrt(group_lengths)*group_norm_sgs[order(sqrt(group_lengths)*group_norm_sgs, decreasing=TRUE)]
  pen_g = sum(gslope_seq*group_norm_sgs) 
  beta = beta[order(abs(beta),decreasing=TRUE)]
  pen_s = sum(slope_seq*abs(beta))
  out = loss + pen_g + pen_s
  return(out)
}

fdr_sensitivity = function(fitted_ids, true_ids,num_coef){ 
  # calculates FDR, FPR, and sensitivity
  num_true = length(intersect(fitted_ids,true_ids))
  num_false = length(fitted_ids) - num_true
  num_missed = length(true_ids) - num_true
  num_true_negatives = num_coef - length(true_ids)
  out=c()
  out$fdr = num_false / (num_true + num_false)
  if (is.nan(out$fdr)){out$fdr = 0}
  out$sensitivity = num_true / length(true_ids)
  if (length(true_ids) == 0){
    out$sensitivity = 1
  }
  out$fpr = num_false / num_true_negatives
  out$f1 = (2*num_true)/(2*num_true + num_false + num_missed)
  if (is.nan(out$f1)){out$f1 = 1}
  return(out)
}

plot_path <- function(beta_matrix, lambdas, how_many,main){
  # Plots the fitted beta values along a lambda path
  beta_order = order(abs(beta_matrix[,dim(beta_matrix)[2]]),decreasing=TRUE)
  max_y = max(beta_matrix)/0.99
  min_y = min(beta_matrix)/0.99
  cols = rainbow(n=how_many)
  plot(x=-log(lambdas), y=beta_matrix[beta_order[1],],type='l',col=cols[1],ylim=c(min_y,max_y),lwd=2,xlab=expression(-log(lambda)),ylab="Fitted value",main=main)

  for (i in 2:how_many){
    lines(x=-log(lambdas),y=beta_matrix[beta_order[i],], type='l', col=cols[i],lwd=2)
  }
}

generate_lambda_path <- function(X,y,groups,alpha,min_frac,nlambda, v_weights, w_weights,group.sizes){
  wts = sort(sqrt(group.sizes),decreasing=TRUE)
  w_weights_expanded = rep(0,length(v_weights))
  group_count = 0
  counter_temp = 1
  for (w_id in order(sqrt(group.sizes),decreasing=TRUE)){
    group_size_temp = group.sizes[w_id]
    w_weights_expanded[(group_count+1):(group_count+group_size_temp)] = rep(w_weights[counter_temp],group_size_temp)*wts[counter_temp]
    group_count = group_count + group_size_temp
    counter_temp = counter_temp+1
  }

  f_0 = mse_grad(y=y,X=X,input=rep(0,dim(X)[2]),num_obs=dim(X)[1])
  f_0_sorted = sort(abs(f_0),decreasing=TRUE)
  
  max_lambda = max(cumsum(f_0_sorted)/cumsum(alpha*v_weights + (1-alpha)*w_weights_expanded)) 

  # calculate the lambda path, including finding lambda_max.

  max_lambda = max_lambda/0.95 # to ensure first variable doesn't enter (first model is null model)
  min_lambda <- min_frac*max_lambda
  
  lambdas = exp(seq(log(max_lambda),log(min_lambda), (log(min_lambda) - log(max_lambda))/(nlambda-1))) 
  #lambdas <- exp(seq(log(max_lambda),log(min_lambda), (log(min_lambda) - log(max_lambda))/(nlambda-1))) # exponential quick decrease sequence (from Simon paper)
  #lambdas = seq(from = max_lambda, to = min_lambda,length.out=nlambda) # linear sequence
  #lambdas = seq(sqrt(max_lambda),sqrt(min_lambda),length.out=nlambda)^2 # quadratic sequence
  #lambdas = sqrt(seq((max_lambda)^2,(min_lambda)^2,length.out=nlambda)) # inverse quadratic

#lambdas <- (seq((max_lambda)^exp(1),(min_lambda)^exp(1), length.out=nlambda))^(1/exp(1)) # inverse exponential sequence

return(lambdas)
}

L1_prox <- function(input, lambda){ # Lasso proximal operator
  out = sign(input) * pmax(0, abs(input) - lambda)
  return(out)
}

group_L1_prox = function(input,lambda,group_info){ # group lasso proximal operator - also has the 0/0 issue
  n_groups = length(unique(group_info))
  out = rep(0,length(input))
  for (i in 1:n_groups){
    grp_idx = which(group_info == unique(group_info)[i])
    if (lambda == 0 & norm(input[grp_idx],type="2") == 0){ # 0/0 = 0
      out[grp_idx] = 0
    } else {
      out[grp_idx] = max((1-(lambda/norm(input[grp_idx],type="2"))),0) * input[grp_idx]}
  }
  return(out)
}

# lambdas of Theorem 2.5 and equation (2.16) in Brzyski et. al. (2016) - from grpSLOPE package
lambdaChiOrtho <- function(fdr, n.group, group.sizes, wt, method) {
  lambda.max <- rep(NA, n.group)
  lambda.min <- rep(NA, n.group)

  for (i in 1:n.group) {
    qchi.seq <- rep(NA, n.group)
    for (j in 1:n.group) {
      qchi.seq[j] <- sqrt(qchisq(1 - fdr*i/n.group, df=group.sizes[j])) / wt[j]
    }
    lambda.max[i] <- max(qchi.seq)
    lambda.min[i] <- min(qchi.seq)
  }

  # stop here if method is "max"
  if (method=="max") return(lambda.max)

  cdfMean <- function(x) {
    pchi.seq <- rep(NA, n.group)
    for (i in 1:n.group) {
      pchi.seq[i] <- pchisq((wt[i]*x)^2, df=group.sizes[i])
    }
    return(mean(pchi.seq))
  }

  lambda.mean <- rep(NA, n.group)
  for (k in 1:n.group) {
    if (lambda.min[k] == lambda.max[k]) {
      lambda.mean[k] <- lambda.max[k]
    } else {
      # compute inverse of cdfMean at 1-fdr*k/n.group
      cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-fdr*k/n.group)),
                             lower = lambda.min[k], upper = lambda.max[k], extendInt="yes")
      lambda.mean[k] <- cdfMean.inv$root
    }
  }

  return(lambda.mean)
}

sgs_var_penalty <- function(q, pen_g,p,lambda,alpha,m,group.sizes,method) {
lambda.max <- rep(0, m)
lambda.min <- rep(0, m)
group.sizes=sort(group.sizes,decreasing=TRUE)
for (i in 1:p){ 
    p_sequence = rep(0,m)
    for (j in 1:m){
        p_sequence[j] = (qnorm(1-(i*q)/(2*p)) - (1/3)*floor(alpha*group.sizes[j])*(1-alpha)*lambda*pen_g[j])/(alpha*lambda)  
    }
    lambda.max[i] <- max(p_sequence)
    lambda.min[i] <- min(p_sequence)
}
if (method == "mean"){
cdfMean <- function(x) {
    p.seq <- rep(0, m)
    for (i in 1:m) {
        p.seq[i] <- pnorm((alpha*lambda*x+(1/3)*floor(alpha*group.sizes[j])*(1-alpha)*lambda*pen_g[i]))
    }
    return(mean(p.seq))
}

lambda.mean <- rep(0, p)
for (k in 1:p) {
    if (lambda.min[k] == lambda.max[k]) {
        lambda.mean[k] <- lambda.max[k]
    } else {
        # compute inverse of cdfMean at 1-fdr*k/n.group
        cdfMean.inv <- uniroot(function(y) (cdfMean(y) - (1-q*k/(2*p))),
                               lower = lambda.min[k], upper = lambda.max[k], extendInt="yes")
        lambda.mean[k] <- max(0,cdfMean.inv$root)
    }
}
return(lambda.mean)
} else if (method == "max"){
  return (lambda.max)
} else {stop("method not valid")}
} 

standardise_sgs <- function(X,y,standardise, intercept,num_obs,type="linear"){
  scale_pen = 1
  standardisation_occured = 0
  y_mean = 0 
  X_center = 0
  X_scale = 1

  if (standardise == "l2") { # l2 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
    X_scale = apply(X,2,function(x) norm(x,type="2"))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:dim(X)[2], function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/sqrt(num_obs)
  } else if (standardise == "sd"){ # sd normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
    X_scale = apply(X, 2, function(x) sqrt(sum(x^2)/num_obs))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:dim(X)[2], function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
  } else if (standardise == "l1"){ # l1 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
    X_scale = colSums(abs(X))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:dim(X)[2], function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/num_obs
  } else { standardisation_occured = 0 } # "none" 
  if (intercept) { # center y and X 
    if (type == "linear"){
    y_mean = mean(y)
    y = y - y_mean}
    if (standardisation_occured == 0){
      X_center = apply(X,2,mean)
      X = sapply(1:dim(X)[2], function(i) X[,i] - X_center[i])
    }
  }

out=c()
out$X=X
out$X_scale = X_scale
out$X_center = X_center
out$y = y
out$y_mean = y_mean
out$scale_pen = scale_pen
return(out)
}

# Compute the usual unbiased estimate of the variance in a linear model. From SLOPE package
estimateNoise <- function(X, y, intercept = TRUE) {
  n <- NROW(X)
  p <- NCOL(X)

  stopifnot(n > p)

  fit <- stats::lm.fit(X, y)
  sqrt(sum(fit$residuals^2) / (n - p + intercept))
}

which_groups <- function(beta, groups){
# outputs the non-zero group ids and effects from beta values
  num_groups = length(unique(groups))
  group.effects = data.frame(group_id = unique(sort(groups)), effect = rep(0,num_groups))
  grp_counter = 1
  for (group_id in unique(groups)){
    group_inds = which(groups==group_id)
    group.effects[grp_counter,]$effect = norm(beta[group_inds], type="2") 
    grp_counter = grp_counter+1 
  }
  selected_group = group.effects[which(group.effects$effect!=0),]$group_id
  group.effects = as(group.effects$effect, "CsparseMatrix")
  rownames(group.effects) = paste0("G", 1:num_groups)
  return(list(selected_group,group.effects))
}

generate_penalties_2 <- function(gFDR, vFDR, pen_method,num_groups,wt_per_grp, len_each_grp, num_vars,alpha,lambda){
  if (pen_method == 1){ # SGS variable mean
    pen_slope_org = BH_sequence(q=vFDR,p=num_vars)
    pen_gslope_org = grp_pen_v1(q=gFDR,pen_v=pen_slope_org,repeats=1e5,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp)
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="max")
  } else if (pen_method == 2){ # SGS variable max
    pen_slope_org = BH_sequence(q=vFDR,p=num_vars)
    pen_gslope_org = grp_pen_adjust(q=gFDR,pen_v=pen_slope_org,repeats=1e5,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp)[[1]]
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="mean")
  } else {stop("method choice not valid")}
out=c()
out$pen_slope_org = pen_slope_org
out$pen_gslope_org = pen_gslope_org
return(out)
}

grp_pen_v1 <- function(q, pen_v,repeats,lambda,alpha,m,group.sizes) {
lambda.max <- rep(0, m)
lambda.min <- rep(0, m)
num_repeats = repeats/length(group.sizes)
z = c()
for (j in 1:num_repeats){
for (i in 1:length(group.sizes)){
    z =append(z,sum(abs(rnorm(group.sizes[i]))))
}
}
num_vars = length(pen_v)
for (i in 1:m){ 
      p_sequence = rep(0,m)
      quantile_v = quantile(z,1-(i*q)/(m))
      for (j in 1:m){
        p_sequence[j] = (quantile_v - alpha*lambda*sum(pen_v[(num_vars - group.sizes[j] + 1):num_vars]))/((1-alpha)*lambda*group.sizes[j])   # pick last few penalties 
    }
    lambda.max[i] <- max(0,max(p_sequence))
}
  return(lambda.max)
} 

run_atos_log_inter <- function(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, step_size, x0, z, u,
                          pen_slope, pen_gslope, f, f_grad,f_opts, f_grad_opts, crossprod_mat,verbose,groups_pen){
  # set values
  wt = as.numeric(sqrt(rep(table(groups),table(groups)))) # Adjust for group weights
  inv_wt = 1/wt
  group_ids = getGroupID(groups) 
  group_ids_pen = getGroupID(groups_pen)
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  num_groups = length(unique(groups))
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  # initial fitting values
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
  prox_input = x0*wt
  z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups_pen,group_id=group_ids_pen)
  z = c(prox_input[1],z)
  z = inv_wt*z
  fz = f(y=y, X=X, input = z, num_obs = num_obs)
  grad_fz = f_grad(y=y, X=X, input = z, num_obs = num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  prox_input = z - (step_size * (grad_fz)) 
  x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
  x = c(prox_input[1],x)
  # fitting
  for (it in 1:max_iter){
    fz = f(y=y, X=X, input = z, num_obs = num_obs)
    grad_fz = f_grad(y=y, X=X, input =z, num_obs = num_obs)
    prox_input = z - (step_size * (u + (grad_fz))) 
    x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
    x = c(prox_input[1],x)
    incr = x - z
    norm_incr = norm(incr,type="2")
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        prox_input = z - (step_size * (u + (grad_fz)))
        x = sortedL1Prox(x=prox_input[-1],lambda=pen_slope*step_size)
        x = c(prox_input[1],x)
        incr = x - z
        norm_incr = norm(incr,type="2")
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y=y, X=X, input = x, num_obs = num_obs) - rhs        
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     
    prox_input = wt*x + (step_size/wt)*u
    z = proxGroupSortedL1(y=prox_input[-1], lambda=pen_gslope*step_size,group=groups_pen,group_id=group_ids_pen)
    z = c(prox_input[1],z)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size

    if (certificate < tol){ # Check for convergence
      success = 1
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = success
out$certificate = certificate
out$it = it
return(out)
}

run_atos <- function(y,X,num_obs, num_vars, groups, backtracking, max_iter, max_iter_backtracking, tol, step_size, x0, z, u,
                          pen_slope, pen_gslope, f, f_grad, f_opts, f_grad_opts, crossprod_mat,verbose){
  # set values
  wt = as.numeric(sqrt(rep(table(groups),table(groups)))) # Adjust for group weights
  inv_wt = 1/wt
  group_ids = getGroupID(groups) 
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  if (is.null(x0)) {x0 = rep(0,num_vars)}
  num_groups = length(unique(groups))
  success = 0 # checks whether convergence happened
  LS_EPS = .Machine$double.eps # R accuracy

  # initial fitting values
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
  z = proxGroupSortedL1(y=x0*wt, lambda=pen_gslope*step_size,group=groups,group_id=group_ids)
  z = inv_wt*z
  fz = f(y=y, X=X, input = z, num_obs = num_obs)
  grad_fz = f_grad(y=y, X=X, input = z, num_obs = num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  x = sortedL1Prox(x=z - (step_size * (grad_fz)) ,lambda=pen_slope*step_size)

  # fitting
  for (it in 1:max_iter){
    fz = f(y=y, X=X, input = z, num_obs = num_obs)
    grad_fz = f_grad(y=y, X=X, input =z, num_obs = num_obs)
    x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))) ,lambda=pen_slope*step_size)
    incr = x - z
    norm_incr = norm(incr,type="2")
    if (norm_incr > 1e-7){
      for (it_ls in 1:max_iter_backtracking){ # Line search
        x = sortedL1Prox(x=z - (step_size * (u + (grad_fz))) ,lambda=pen_slope*step_size)
        incr = x - z
        norm_incr = norm(incr,type="2")
        rhs = fz + crossprod_mat(grad_fz,incr) + (norm_incr ^ 2) / (2 * step_size)
        ls_tol = f(y=y, X=X, input = x, num_obs = num_obs) - rhs        
        if (as.numeric(ls_tol) <= as.numeric(LS_EPS)){
          break
        }
        else {
          step_size = step_size*backtracking 
        }
      }
    }     

    z = proxGroupSortedL1(y=wt*x + (step_size/wt)*u, lambda=pen_gslope*step_size,group=groups,group_id=group_ids)
    z = z*inv_wt
    u = u + (x - z) / step_size

    certificate = norm_incr / step_size

    if (certificate < tol){ # Check for convergence
      success = 1
      break
    }
  if (verbose==TRUE){print(paste0("Iteration: ", it,"/",max_iter, " done"))}
  }
out=c()
out$x = x
out$u = u
out$z = z
out$success = success
out$certificate = certificate
out$it = it
return(out)
}