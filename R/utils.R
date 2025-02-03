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
# Stop 1 dim matrices becoming vectors
old <- `[`
`[` <- function(...) { old(..., drop=FALSE) }
  
# -------------------------------------------------------------
# generalf functions
# -------------------------------------------------------------
# checks if a variable is binary
is.binary <- function(v) {
    x <- unique(v)
    length(x) - sum(is.na(x)) == 2L
}

is.decreasing <- function(x){ # checks if a sequence is decreasing
  if (length(x) == 1){
    state = TRUE
  } else {
  state = TRUE
  for (i in 2:length(x)){
    if (x[i] > x[i-1]){
      state = FALSE
    }
  }
  }
  return(state)
}

is.strictly.decreasing <- function(x){ # checks if a sequence is decreasing
  if (length(x) == 1){
    state = TRUE
  } else {
  state = TRUE
  for (i in 2:length(x)){
    if (x[i] >= x[i-1]){
      state = FALSE
    }
  }
  }
  return(state)
}

sigmoid = function(x) {
   1 / (1 + exp(-x))
}

getGroupID <- function(group) { # from grpSLOPE package, which is no longer available on CRAN
  group.unique <- unique(group)
  n.group <- length(group.unique)
  group.id <- list()
  for (i in 1:n.group){
    id <- as.character(group.unique[i])
    group.id[[id]] <- which(group==group.unique[i])
  }
  return(group.id)
}

norm_vec <- function(x) sqrt(sum(x^2))

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
    group_norm_sgs[i] = norm_vec(beta[(len_so_far+1):(len_so_far+g_length)])
    len_so_far = len_so_far+g_length
  }
  if (intercept==TRUE){
    loss = mse_loss(y=y,Xbeta=arma_mv(Matrix::cbind2(1,X),beta_org),num_obs = num_obs, crossprod_mat = base::crossprod)
  } else {loss = mse_loss(y=y,Xbeta=arma_mv(X,beta_org),num_obs = num_obs, crossprod_mat = base::crossprod)}
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

# lambdas of Theorem 2.5 and equation (2.16) in Brzyski et al. (2016) - from grpSLOPE package, which is no longer available on CRAN
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

standardise_data <- function(X,y,standardise, intercept,num_obs,type="linear"){
  scale_pen = 1
  standardisation_occured = 0
  y_mean = 0
  X_center = rep(0,ncol(X))
  X_scale = rep(1,ncol(X))

  if (standardise == "l2") { # l2 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = apply(X,2,function(x) norm_vec(x))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/sqrt(num_obs)
  } else if (standardise == "sd"){ # sd normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = apply(X, 2, function(x) sqrt(sum(x^2)/num_obs))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
  } else if (standardise == "l1"){ # l1 normalisation
    X_center = apply(X,2,mean)
    X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
    X_scale = colSums(abs(X))
    if (any(X_scale==0)){
      stop("not able to standardise X as there exists at least one predictor with no variance")
    }
    X = sapply(1:ncol(X), function(i) X[,i]/X_scale[i])
    standardisation_occured = 1
    scale_pen = 1/num_obs
  } else { standardisation_occured = 0 } # "none"
  if (intercept) { # center y and X
    if (type == "linear"){
    y_mean = mean(y)
    y = y - y_mean}
    if (standardisation_occured == 0 & !inherits(X,"dgCMatrix")){
      X_center = apply(X,2,mean)
      X = sapply(1:ncol(X), function(i) X[,i] - X_center[i])
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

which_groups <- function(beta, groups){
# outputs the non-zero group ids and effects from beta values
  num_groups = length(unique(groups))
  group_effects = data.frame(group_id = unique(sort(groups)), effect = rep(0,num_groups))
  grp_counter = 1
  for (group_id in unique(groups)){
    group_inds = which(groups==group_id)
    group_effects[grp_counter,]$effect = norm_vec(beta[group_inds])
    grp_counter = grp_counter+1
  }
  selected_grp = group_effects[which(group_effects$effect!=0),]$group_id
  group_effects = as.vector(group_effects$effect)
  return(list(selected_grp,group_effects))
}

l2_group_operator = function(x,P, groupIDs,power){
    out = rep(0,length(groupIDs))
    for (g in 1:length(groupIDs)){
        out[g] = (P[g]^power)*norm_vec(x[unlist(groupIDs[g])])
    }
    return(out)
}

soft_thresholding_operator <- function(x,thres){
    out = sign(x)*ifelse(abs(x) - thres <=0,0, abs(x) - thres)
    return(out)
}

path_shape <- function(lambda_max, path_length, min_frac){
  min_lambda = min_frac*lambda_max
  lambda_seq = exp(seq(log(lambda_max),log(min_lambda), (log(min_lambda) - log(lambda_max))/(path_length-1))) 
  return(lambda_seq)
}

reorder_group <- function(groups){
  max_grp_id = length(unique(groups))
  new_grp = rep(0,length(groups))
  all_grp_indices = as.numeric(names(table(groups)))
  for (i in 1:max_grp_id){
	  var_ids = which(groups == all_grp_indices[i]) 
	  new_grp[var_ids] = i
  }
return(new_grp)
}

# -------------------------------------------------------------
# loss and grad functions
# -------------------------------------------------------------
### gaussian
mse_loss <- function(y, Xbeta, num_obs, crossprod_mat){ # linear loss function
  return(as.double(crossprod_mat(y-Xbeta)/(2*num_obs)))
}

mse_grad <- function(y, Xbeta, num_obs){ # pseudo-gradient of loss, need to multiply by X^T
    return((Xbeta-y)/num_obs)
}

### binomial
log_loss <- function(y,Xbeta,num_obs, crossprod_mat){ # logistic loss function. stable version for y{0,1}
  return(mean((1-y)*Xbeta - logsig(Xbeta)))
}

log_grad <- function(y,Xbeta,num_obs){# stable version for y{0,1}. pseudo-gradient of loss, need to multiply by X^T
  return(expit_b(Xbeta,y)/num_obs)
}

# stable log implementations - from https://fa.bianp.net/blog/2019/evaluate_logistic/
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

# -------------------------------------------------------------
# path functions
# -------------------------------------------------------------
gen_path_sgs = function(grad_vec_zero, groups, groupIDs, alpha, w_weights, v_weights, path_length, min_frac, group_sizes=NULL){
  num_vars = length(groups)
  grad_sorting = order(abs(grad_vec_zero), decreasing=TRUE)
  grad_sorting_groups = l2_group_operator(x=grad_vec_zero,P=table(groups),groupIDs,-0.5)
  group_sorting = order(grad_sorting_groups,decreasing=TRUE)
  top = rep(0,num_vars)
  bottom = rep(0,num_vars)
  for (g in 1:num_vars){ # try also ordering by gradient
      grp_index = match(groups[g],group_sorting)
      size_pen = sqrt(table(groups))[groups[g]]
      var_index = match(g,grad_sorting)
      top[g] = abs(grad_vec_zero[g])
      bottom[g] = alpha*v_weights[var_index] + size_pen*(1-alpha)*w_weights[grp_index]
  }
  lambda_max = max(cumsum(sort(top,decreasing=TRUE))/cumsum(sort(bottom,decreasing=TRUE)))/0.99
  lambda_seq = path_shape(lambda_max,path_length,min_frac)
  return(lambda_seq)
}

gen_path_gslope = function(grad_vec_zero, groups, groupIDs, alpha=0, w_weights, v_weights=NULL, path_length, min_frac, group_sizes=NULL){
  l2_vector = l2_group_operator(x=grad_vec_zero, P=table(groups), groupIDs, power=-0.5)
  lambda_max = max(cumsum(sort(abs(l2_vector),decreasing=TRUE))/cumsum(w_weights))/0.99
  lambda_seq = path_shape(lambda_max,path_length,min_frac)
  return(lambda_seq)
}

# -------------------------------------------------------------
# penalty functions
# -------------------------------------------------------------
BH_sequence = function(q,p){
  # Calculates the BH sequences for SLOPE and gSLOPE
  p_sequence = rep(0,p)
  for (i in 1:p){
    q_i = (i*q)/(2*p)
    p_sequence[i] = qnorm(1-q_i)
  }
return(p_sequence)
}

sgs_var_penalty <- function(q, pen_g, p, lambda, alpha, m, group.sizes, method){
    lambda.max = rep(0, m)
    lambda.min = rep(0, m)
    group.sizes = sort(group.sizes, decreasing = TRUE)
    for (i in 1:p) {
        p_sequence = rep(0, m)
        for (j in 1:m) {
            p_sequence[j] = (qnorm(1 - (i * q) / (2 * p)) - (1 / 3) * floor(alpha * group.sizes[j]) * (1 - alpha) * lambda * pen_g[j]) / (alpha * lambda)
        }
        lambda.max[i] <- max(p_sequence)
        lambda.min[i] <- min(p_sequence)
    }
    if (method == "mean") {
        cdfMean <- function(x) {
            p.seq <- rep(0, m)
            for (i in 1:m) {
                p.seq[i] <- pnorm((alpha * lambda * x + (1 / 3) * floor(alpha * group.sizes[j]) * (1 - alpha) * lambda * pen_g[i]))
            }
            return(mean(p.seq))
        }
        
        lambda.mean <- rep(0, p)
        for (k in 1:p) {
            if (lambda.min[k] == lambda.max[k]) {
                lambda.mean[k] = lambda.max[k]
            } else {
                # compute inverse of cdfMean at 1-fdr*k/n.group
                cdfMean.inv = uniroot(function(y) (cdfMean(y) - (1 - q * k / (2 * p))), lower = lambda.min[k], upper = lambda.max[k], extendInt = "yes")
                lambda.mean[k] = max(0, cdfMean.inv$root)
            }
        }
        return(lambda.mean)
    } else if (method == "max") {
        return(lambda.max)
    } else {
        stop("method not valid")
    }
}

### as-sgs penalty functions
gen_pens_as_sgs <- function(gFDR, vFDR, pen_method,groups,alpha,lambda){
  num_vars = length(groups)
  group_ids = getGroupID(groups)
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  wt_per_grp = wt_per_grp[names(group_ids)]
  num_groups = length(unique(groups))
  if (pen_method == 1){ # SGS variable mean
    pen_slope_org = BH_sequence(q=vFDR,p=num_vars)
    pen_gslope_org = grp_pen_as_sgs_mean(q=gFDR,pen_v=pen_slope_org,repeats=1e5,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp)
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="max")
  } else if (pen_method == 2){ # SGS variable max
    pen_slope_org = BH_sequence(q=vFDR,p=num_vars)
    pen_gslope_org = grp_pen_as_sgs_max(q=gFDR,pen_v=pen_slope_org,repeats=1e5,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp)[[1]]
    pen_slope_org = sgs_var_penalty(q=vFDR, pen_g=pen_gslope_org,p=num_vars,lambda=lambda,alpha=alpha,m=num_groups,group.sizes=len_each_grp,method="mean")
  } else {stop("method choice not valid")}
  out=c()
  out$pen_slope_org = pen_slope_org
  out$pen_gslope_org = pen_gslope_org
  return(out)
}

grp_pen_as_sgs_mean <- function(q, pen_v, repeats, lambda, alpha, m, group.sizes){
  lambda.max = rep(0, m)
  lambda.min = rep(0, m)
  num_repeats = repeats / length(group.sizes)
  z = c()
  for (j in 1:num_repeats) {
    for (i in 1:length(group.sizes)) {
      z = append(z, sum(abs(rnorm(group.sizes[i]))))
    }
  }
  num_vars = length(pen_v)
  for (i in 1:m) {
    p_sequence = rep(0, m)
    quantile_v = quantile(z, 1 - (i * q) / (m))
    for (j in 1:m) {
      p_sequence[j] = (quantile_v - alpha * lambda * sum(pen_v[(num_vars - group.sizes[j] + 1):num_vars])) / ((1 - alpha) * lambda * group.sizes[j])   # pick last few penalties
    }
    lambda.max[i] = max(0, max(p_sequence))
  }
  return(lambda.max)
}

grp_pen_as_sgs_max <- function(q, pen_v, repeats, lambda, alpha, m, group.sizes){
    lambda.max = rep(0, m)
    lambda.min = rep(0, m)
    num_repeats = repeats / length(group.sizes)
    z = c()
    for (j in 1:num_repeats) {
        for (i in 1:length(group.sizes)) {
            z = append(z, sum(abs(rnorm(group.sizes[i]))))
        }
    }
    group.sizes = sort(group.sizes, decreasing = TRUE)
    num_vars = length(pen_v)
    z_fn = ecdf(z)
    an.error.occured = FALSE
    successful_roots = 0
    attempts = 0
    counter = 0
    counter2 = 0
    while (successful_roots != m & attempts < 10) {
        an.error.occured = FALSE
        successful_roots = 0
        for (i in 1:m) {
            counter = 0
            p_sequence = rep(0, m)
            quantile_v = quantile(z, 1 - (i * q) / (m))
            for (j in 1:m) {
                p_sequence[j] = (quantile_v - alpha * lambda * sum(pen_v[(counter + 1):(counter + group.sizes[j])])) / ((1 - alpha) * lambda * group.sizes[j])
                counter = counter + group.sizes[j]
            }
            lambda.max[i] = max(p_sequence)
            lambda.min[i] = min(p_sequence)
        }
        cdfMean <- function(x) {
            p.seq = rep(0, m)
            counter2 = 0
            for (i in 1:m) {
                p.seq[i] = z_fn(sum(abs(x)) * (1 - alpha) * lambda * group.sizes[j] + alpha * lambda * sum(pen_v[(counter2 + 1):(counter2 + group.sizes[j])]))
                counter2 = counter2 + group.sizes[j]
            }
            return(mean(p.seq))
        }
        
        lambda.mean <- rep(0, m)
        for (k in 1:m) {
            if (lambda.min[k] == lambda.max[k]) {
                lambda.mean[k] = lambda.max[k]
            } else {
                # compute inverse of cdfMean at 1-fdr*k/n.group
                tryCatch({
                    cdfMean.inv = uniroot(function(y) (cdfMean(y) - (1 - q * k / (m))), lower = lambda.min[k], upper = lambda.max[k], extendInt = "yes")
                }, error = function(e) {
                    an.error.occured <<- TRUE
                })
                if (an.error.occured) {
                    
                } else {
                    successful_roots = successful_roots + 1
                    lambda.mean[k] = max(0, cdfMean.inv$root)
                }
            }
        }
        if (an.error.occured) {
            pen_v = pen_v * 0.95
            attempts = attempts + 1
        }
    }
    if (attempts >= 9) {
        print("error")
    }
    return(list(lambda.mean, pen_v))
}

# -------------------------------------------------------------
# screen functions
# -------------------------------------------------------------
inner_screening = function(c, lambda, lambda_new){ # algorithm 2
  m = length(lambda)
  c_input = c
  c = sort(abs(c), decreasing=TRUE) + lambda - 2*lambda_new 
  i = 1
  k = 0
  s = 0
  while (i+k <= m){
      s = s + c[i+k]
      if (s >= 0){
          k = k+i
          i = 1
          s = 0
      } else {
          i = i+1
      }
  }
  output = order(abs(c_input),decreasing=TRUE)[1:k] # corresponds to S
  return(output)
}

## sgs
# screening rules
sgs_grp_screen <- function(grad_vec, current_beta, tbl_grps, groupIDs, alpha, pen_slope_org, pen_gslope_org, lambda_new, lambda, wt=NULL){
  beta_order = order(abs(current_beta))
  soft_input = soft_thresholding_operator(x=grad_vec[beta_order],thres=lambda*pen_slope_org*alpha)
  soft_input[beta_order] = soft_input
  input_screen = l2_group_operator(x=soft_input,P=tbl_grps, groupIDs, power=-0.5) # do we sort grad_vec? I think so
  screen_set_grp = inner_screening(input_screen, lambda=lambda*pen_gslope_org*(1-alpha), lambda_new = lambda_new*pen_gslope_org*(1-alpha)) # corresponds to S  
  return(screen_set_grp)
}

sgs_var_screen <- function(grad_vec, groupIDs, screen_set_grp, alpha, pen_slope_org, lambda_new, lambda, active_set_var=NULL){
  screen_set_var_initial = unlist(groupIDs[screen_set_grp])
  screen_set_var = inner_screening(grad_vec[screen_set_var_initial], lambda=alpha*lambda*pen_slope_org[1:length(screen_set_var_initial)], lambda_new = alpha*lambda_new*pen_slope_org[1:length(screen_set_var_initial)]) # corresponds to S
  screen_set_var = screen_set_var_initial[screen_set_var]
  return(screen_set_var)
}

# kkt check
sgs_kkt_check = function(grad_vec, current_beta, groups, groupIDs, alpha, pen_slope_org, pen_gslope_org, lambda, tbl_grps, machine_tol, epsilon_set_var=NULL, non_zero_groups){ # selects only those penalties associated with the zero groups, so the bottom ones - think this is the correct approach
  # group check
  if (length(non_zero_groups) == 0){
    zero_groups = unique(groups)
  } else {
    zero_groups = which(!(unique(groups) %in% non_zero_groups))
  }
    
  beta_order = order(abs(current_beta),decreasing=TRUE)
  soft_input = soft_thresholding_operator(x=grad_vec[beta_order],thres=lambda*pen_slope_org*alpha)
  soft_input[beta_order] = soft_input 
  inner_val = l2_group_operator(x=soft_input,P=tbl_grps, groupIDs, power=-0.5) # check whether pen_seq_new or pen_seq?
  sort_index_gbn = order(inner_val,decreasing=TRUE)[1:length(which(inner_val!=0))]
  sort_index_gbn = sort_index_gbn[sort_index_gbn %in% zero_groups]
  
  subdiff = cumsum(sort(inner_val[sort_index_gbn],decreasing=TRUE) - (1-alpha)*lambda*pen_gslope_org[1:length(sort_index_gbn)])
  violations_grp = rep(0,length(tbl_grps))
  violation_ids = which(subdiff > machine_tol)
  violations_grp[sort_index_gbn[violation_ids]] = 1
  violations_grp[non_zero_groups] = 0  

  # var check - for variables within groups that are not zero (including violations)
  grps_to_check = union(non_zero_groups,which(violations_grp==1))
  vars_to_check = unlist(groupIDs[grps_to_check])

  sort_index_var = order(abs(grad_vec[vars_to_check]),decreasing=TRUE)
  subdiff = cumsum(abs(grad_vec[vars_to_check][sort_index_var]) - alpha*lambda*pen_slope_org[1:length(vars_to_check)])
  
  violations_var = rep(0,length(current_beta))
  violation_ids = which(subdiff > machine_tol)
  violations_var[vars_to_check[sort_index_var[violation_ids]]] = 1

  return(list(violations_var,violations_grp))
}

## gslope
# screening rule
gslope_grp_screen <- function(grad_vec, current_beta=NULL, tbl_grps, groupIDs, alpha=0, pen_slope_org=NULL, pen_gslope_org, lambda_new, lambda, wt= NULL){
  screen_set_grp = inner_screening(l2_group_operator(x=grad_vec, P=tbl_grps, groupIDs, power=-0.5), lambda=lambda*pen_gslope_org, lambda_new = lambda_new*pen_gslope_org) # corresponds to S
  return(screen_set_grp)
}

# kkt check
gslope_kkt_check = function(grad_vec, current_beta, groups, groupIDs, alpha=0, pen_slope_org=NULL, pen_gslope_org, lambda, tbl_grps, machine_tol, epsilon_set_var=NULL, non_zero_groups){ # selects only those penalties associated with the zero groups, so the bottom ones - think this is the correct approach
  grad_group_operator = l2_group_operator(x=grad_vec,P=tbl_grps,power=-0.5,groupIDs=groupIDs)
  if (length(non_zero_groups) == 0){
    zero_groups = unique(groups)
  } else {
    zero_groups = which(!(unique(groups) %in% non_zero_groups))
  }
  
  sort_index_gbn = order(grad_group_operator,decreasing=TRUE)
  sort_index_gbn = sort_index_gbn[sort_index_gbn %in% zero_groups]
  
  subdiff = cumsum(sort(grad_group_operator[sort_index_gbn],decreasing=TRUE) - lambda*pen_gslope_org[(length(non_zero_groups)+1):length(tbl_grps)])

  violations = rep(0,length(tbl_grps))
  violation_ids = which(subdiff > machine_tol)
  violations[sort_index_gbn[violation_ids]] = 1
  violations[non_zero_groups] = 0  
return(list(1,violations))
}

# -------------------------------------------------------------
# algorithm functions
# -------------------------------------------------------------
init_lipschitz <- function(f, f_grad, mult_fcn, x0, X, y, num_obs, tX, crossprod_mat){
  L0 = 1e-3
  Xx0 = mult_fcn(X,x0)
  f0 = f(y, Xx0, num_obs, crossprod_mat)
  grad0 = mult_fcn(tX,f_grad(y, Xx0, num_obs))

  x_tilde = x0 - (1 / L0)*grad0
  f_tilde = f(y, mult_fcn(X,x_tilde), num_obs, crossprod_mat) 

  for (i in 1:100){
    if (f_tilde <= f0){
      break
    } else {
      L0 = L0 * 10
      x_tilde = x0 - (1 / L0) * grad0
      f_tilde = f(y, mult_fcn(X,x_tilde), num_obs, crossprod_mat) 
    }
  }
  return(L0)
}

proxGroupSortedL1 <- function(y, lambda,group, group_id, num_groups) {
  # proximal operator for group SLOPE - adapted so that the 0/0 = NaN error doesn't occur
  # adapted from grpSLOPE package, which is no longer available on CRAN
  if (length(lambda) != num_groups) {
    stop("Length of lambda should be equal to the number of groups.")
  }

  # compute Euclidean norms for groups in y
  group_norm <- rep(NA, num_groups)
  for (i in 1:num_groups){
    selected <- group_id[[i]]
    group_norm[i] <- norm_vec(y[selected])
  }

  # get Euclidean norms of the solution vector
  prox_norm <- sortedL1Prox(x=group_norm, lambda=lambda, method="stack")

  # compute the solution
  prox_solution <- rep(NA, length(y))
  for (i in 1:num_groups){
    selected <- group_id[[i]]
    if (group_norm[i] == 0){ # to stop 0/0 = NaN
      prox_solution[selected] <- 0
      } else {
    prox_solution[selected] <- prox_norm[i] / group_norm[i] * y[selected] }
  }
  return(prox_solution)
}

# -------------------------------------------------------------
# method functions
# -------------------------------------------------------------
# Compute the usual unbiased estimate of the variance in a linear model. From SLOPE package
estimateNoise <- function(X, y, intercept = TRUE) {
  n <- nrow(X)
  p <- ncol(X)

  stopifnot(n > p)

  fit <- stats::lm.fit(X, y)
  sqrt(sum(fit$residuals^2) / (n - p + intercept))
}

check_group_vector <- function(vec) {
  # Check if the vector is sorted and has no gaps
  is_sorted <- all(diff(vec) >= 0)
  has_no_gaps <- all(diff(unique(vec)) == 1)
  
  return(is_sorted && has_no_gaps)
}