source("utils.R")
#source("toy_data_gen.R")

fit_sgs <- function(X, y, groups, pen_method = 1, type="linear", lambda, alpha=0.95, vFDR=0.1, gFDR=0.1, max_iter = 5000, backtracking = 0.7, max_iter_backtracking = 100, tol = 1e-5, standardise = "l2", intercept = TRUE,w_weights=NULL,v_weights=NULL,x0 = NULL, u = NULL,verbose=FALSE){
  
  # -------------------------------------------------------------
  # checks
  # ------------------------------------------------------------- 
  if (anyNA(y) | anyNA(X)) {
    stop("input contains missing values")
  }
  if (!(is.matrix(X) | is(X, 'sparseMatrix')) | !(is.matrix(y) | is.vector(y))){
    stop("X and y must be matrices/vectors. Use the as.matrix function to convert")
  }
  if (length(y) == 0) {
    stop("y is empty")
  }
  if (nrow(X) == 0) {
    stop("X is empty")
  }
  if (length(y) != nrow(X)) {
    stop("the number of samples in y must match the number of rows in X")
  }
  if (lambda<0){
    stop("lambda can not be negative")
  }
  if (alpha<0 | alpha>1){
    stop("alpha must be in [0,1]")
  }
  if (vFDR<=0 | vFDR>=1 | gFDR<=0 | gFDR>=1){
    stop("FDR must be in (0,1)")
  }
  if (alpha == 1 & lambda != 0){
    warning("this package is not optimised for SLOPE. Consider using the SLOPE package instead")
  }
  if (alpha == 0 & lambda != 0){
    warning("this package is not optimised for group SLOPE. Consider using the grpSLOPE package instead")
  }
  if (lambda == 0){
    warning("this package is not optimised for OLS. Consider using the lm function instead")
  }

  # -------------------------------------------------------------
  # pre-process data
  # ------------------------------------------------------------- 
  num_vars = dim(X)[2]
  num_obs = dim(X)[1]

  if (sum(X==0) > (num_vars*num_obs)/2){
    warnings("X appears to be a sparse matrix. Try converting to dgCMatrix type for improved performance")
  }
  
  is_sparse = FALSE
  crossprod_mat = base::crossprod
  if (inherits(X,"dgCMatrix")){ # check if matrix is sparse
    crossprod_mat = Matrix::crossprod
    is_sparse = TRUE
  }

  if (standardise!="none" & is_sparse){
    stop("standardising a matrix that is sparse. this would remove sparsity of X. set standardise to none")
  }

  #if (type == "logistic"){ # change to [-1,1] to aid training
  #  y = ifelse(y==1,1,-1)
  #}

  # standardise
  if (standardise=="none"){
      scale_pen = 1
      y_mean = 0 
      X_center = 0
      X_scale = 1
    } else {
      standardise_out = standardise_sgs(X=X,y=y,standardise=standardise,intercept=intercept,num_obs=num_obs,type=type)
      X = standardise_out$X
      X_scale = standardise_out$X_scale
      X_center = standardise_out$X_center
      y = standardise_out$y
      y_mean = standardise_out$y_mean
      scale_pen = standardise_out$scale_pen
    }

  # -------------------------------------------------------------
  # set values
  # ------------------------------------------------------------- 
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

  # type of model
  if (type == "linear"){ 
    if (is_sparse){
      f = mse_loss_sparse
      f_grad = mse_grad_sparse
    } else {
      f = mse_loss
      f_grad = mse_grad
    }
  } else if (type == "logistic"){
    if (is_sparse){
      f = log_loss_sparse
      f_grad = log_grad_sparse
    } else {
      f = log_loss
      f_grad = log_grad
    }
  } else {stop("loss function not supported")} 
  f_opts = list(y=y,X=X,num_obs=num_obs)
  f_grad_opts = list(y=y,X=X,num_obs=num_obs)
  lambda_org = lambda
  lambda = scale_pen*lambda

  # weights
  if (is.null(v_weights) & is.null(w_weights)){
    pens_out = generate_penalties(gFDR, vFDR, pen_method,num_groups,wt_per_grp, len_each_grp, num_vars,alpha)
    pen_slope_org = pens_out$pen_slope_org
    pen_gslope_org = pens_out$pen_gslope_org
    pen_slope = alpha*lambda*pen_slope_org
    pen_gslope = (1-alpha)*lambda*pen_gslope_org
  } else {
    pen_slope_org = v_weights
    pen_gslope_org = w_weights
    pen_slope = alpha*lambda*v_weights
    pen_gslope = (1-alpha)*lambda*w_weights
  } 

  # penalty checks
  if (any(pen_slope < 0) | any(pen_gslope < 0)){
    stop("penalty sequences must be positive")
  }

  if (!is.decreasing(pen_slope) | !is.decreasing(pen_slope)){
    stop("penalty sequences must be decreasing")
  }

  # -------------------------------------------------------------
  # initial fitting values
  # ------------------------------------------------------------- 
  step_size = 1/init_lipschitz(f=f,f_grad=f_grad, x0=x0, f_opts = f_opts, f_grad_opts = f_grad_opts)
 # rand.mat <- matrix(rnorm(num_vars), c(num_vars, 1))
 # rand.mat <- rand.mat / norm(rand.mat, "f")
 # rand.mat <- t(X) %*% (X %*% rand.mat)
 # L <- (1/num_obs)*norm(rand.mat, "f")*wt[1]*wt[1]
 # step_size = 1/L
  
  z = proxGroupSortedL1(y=x0*wt, lambda=pen_gslope*step_size,group=groups,group_id=group_ids)
  z = inv_wt*z
  fz = f(y=y, X=X, input = z, num_obs = num_obs)
  grad_fz = f_grad(y=y, X=X, input = z, num_obs = num_obs) # loss gradient at z
  if (is.null(u)) {u= rep(0,num_vars)}
  x = sortedL1Prox(x=z - (step_size * (grad_fz)) ,lambda=pen_slope*step_size)

  # -------------------------------------------------------------
  # fitting
  # ------------------------------------------------------------- 
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

  if (success == 0){ # check for convergence
    warning("algorithm did not converge, try using more iterations")
  } else {if (verbose==TRUE){print("Algorithm converged")}}

  # -------------------------------------------------------------
  # generate output
  # ------------------------------------------------------------- 
  out = c()
  out$x_beta = x
  if (max((x-z)^2) < 1e-3 & mean((x-z)^2) < 1e-3){ # if solutions are very similar, pick more stable version
    if (length(which(x!=0)) <= length(which(z!=0))){ # Picking the solution with less residual values, if this is true, x is picked
      out$beta = as.matrix(x)
    } else {
      out$beta = as.matrix(z)
    }
  } else { # if solutions aren't similar, pick x
    out$beta = as.matrix(x)
  }

  # scale beta depending on transformations
  if (standardise!="none"){ 
    out$beta = out$beta/X_scale
    out$x_beta = out$x_beta/X_scale
  }

  if (length(out$beta[out$beta!=0]) != 0){
    if (min(abs(out$beta[out$beta!=0])) < 1e-3){ # remove small resid values if solution not stable
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>(1e-4))]))[4]*1e-3
      threshold_x = quantile(abs(out$beta[which(abs(out$beta)>=threshold_x)]))[4]*1e-2
      if (!is.na(threshold_x) & lambda_org!=0) { # When lambda = 0, we don't want to remove small values, as no penalisation is occuring
        threshold_x = ifelse(threshold_x>1e-2,1e-2, threshold_x) # if threshold too big, set to 1e-2 
        out$beta = ifelse(abs(out$beta)>threshold_x,out$beta,0)
        out$beta = as.matrix(out$beta)
      }
    }
  }

  out$selected_var = which(out$beta!=0)
  which_groups_out = which_groups(beta=out$beta,groups=groups)
  out$selected_group = which_groups_out[[1]]
  out$group.effects = which_groups_out[[2]]

  if (intercept){ # get beta back to original scale
    out$beta = as.matrix(c(y_mean - sum(X_center*out$beta),out$beta))
    out$x_beta = as.matrix(c(y_mean - sum(X_center*out$x_beta),out$x_beta))
  } 

  if (is.null(colnames(X))){ # Add variable names to output
    if (intercept){
        rownames(out$beta) = c("(Intercept)", paste0("v", 1:(num_vars)))
      } else {
        rownames(out$beta) = paste0("v", 1:num_vars)
      }
    } else {
    if (intercept){
        rownames(out$beta) = c("(Intercept)", colnames(X))
      } else {
        rownames(out$beta) = colnames(X)
      }
  }
  out$beta = as(out$beta,"CsparseMatrix")
  out$z = z
  out$x = x
  out$u = u
  out$type = type
  out$pen_slope = pen_slope_org
  out$pen_gslope = pen_gslope_org
  out$lambda=lambda_org
  out$success = success
  out$num_it = it
  out$certificate = certificate
  out$intercept = intercept
  class(out) <- "SGS"
  
  return(out)
}

print.SGS = function(out, digits = max(3, getOption("digits") - 3), ...){ 
  num.nonzero <- apply(out$beta,2, function(z){sum(z != 0)})
  cat("\n regression type: ", out$type, "\n\n")
  print(cbind(lambdas = out$lambdas, num.nonzero = num.nonzero, convergence = out$success))
}

predict.SGS = function(fit, predict_X, type="linear"){
  if (type=="linear"){
    if (fit$intercept){
    predictions = arma_mm(cbind(1,predict_X),as.vector(fit$beta))
    } else {
      predictions = arma_mm(predict_X,as.vector(fit$beta))
    }
  } else if (type == "logistic"){
    predictions = c()
    if (fit$intercept){
      predictions$response = sigmoid(arma_mm(cbind(1,predict_X),as.vector(fit$beta)))
    } else {
      predictions$response = sigmoid(arma_mm(predict_X,as.vector(fit$beta)))
    }
    predictions$class = ifelse(predictions$response>0.5,1,0)
  } else {
       stop("not a valid type")
  }
  return(predictions)
}