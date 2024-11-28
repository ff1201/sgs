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

# File contains general wrapper code for fitting models, preparing the inputs and outputs accordingly

general_fit <- function(X, y, groups, model, path_fcn, var_screen_fcn, grp_screen_fcn, kkt_fcn, type, lambda, path_length, alpha, vFDR, gFDR, pen_method, backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, v_weights, w_weights, screen, verbose, gamma_1, gamma_2){
 
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
  if (type == "logistic" & !is.binary(y)){
    stop("when using logistic regression the response, y, needs to be a binary variable")
  }
  if (type == "linear" & is.binary(y)){
    stop("using a binary variable with linear regression. Use logistic instead")
  }
  if (any(lambda<0)){
    stop("lambda can not be negative")
  }
  if (alpha < 0 | alpha > 1){
    stop("alpha must be in [0,1]")
  }
  if (gFDR<=0 | gFDR>=1){
    stop("FDR must be in (0,1)")
  }
  if (type == "logistic" & intercept == TRUE){
    warning("logistic regression with intercept has a bug. See sgs GitHub page for more details")
  }
  if (!check_group_vector(groups)){ # check if group indexes ordered correctly
    warning("group indices are not set from 1, contain gaps, or are not ordered, which can lead to a bug during fitting. See GitHub page for more details.")
  } 

  # identify fit type
  if (any(lambda == "path") | length(lambda) > 1){
    fit_type = "path"
  } else {
    fit_type = "single"
  }

  # -------------------------------------------------------------
  # pre-process data
  # -------------------------------------------------------------
  if (inherits(X,"dgCMatrix")){ # check if matrix is sparse
    crossprod_mat = Matrix::crossprod
    mult_fcn = arma_sparse
    if (standardise!="none"){
      stop("standardising a matrix that is sparse. this would remove sparsity of X. set standardise to none")
    }
  } else {
    crossprod_mat = base::crossprod
    mult_fcn = arma_mv
  }

  # type of model
  if (type == "linear"){ 
    f_grad = mse_grad
    f = mse_loss
  } else if (type == "logistic"){
    f_grad = log_grad
    f = log_loss
  } else {stop("loss function not supported")} 

  groupIDs = getGroupID(groups)
  num_groups = length(groupIDs)
  num_vars = ncol(X)
  num_obs = nrow(X)
  wt = as.numeric(sqrt(rep(table(groups),table(groups))))
  if (type == "logistic" & intercept){num_vars = ncol(X)+1} # fit intercept directly for logistic model

  if (sum(X==0) > (num_vars*num_obs)/2){
    warnings("X appears to be a sparse matrix. Try converting to dgCMatrix type for improved performance")
  }

  # standardise
  if (intercept){ # calculate zero gradient
        grad_vec_zero = mult_fcn(Matrix::t(Matrix::cbind2(1,X)), f_grad(y, mult_fcn(Matrix::cbind2(1,X),rep(0,num_vars+1)), num_obs = num_obs))[-1]
    } else {grad_vec_zero = mult_fcn(Matrix::t(X),f_grad(y, mult_fcn(X,rep(0,num_vars)), num_obs = num_obs))}

  if (standardise=="none" & intercept==FALSE){
      scale_pen = 1
      y_mean = 0 
      X_center = rep(0,num_vars)
      X_scale = rep(1,num_vars)
    } else {
      standardise_out = standardise_data(X=X,y=y,standardise=standardise,intercept=intercept,num_obs=num_obs,type=type)
      X = standardise_out$X
      X_scale = standardise_out$X_scale
      X_center = standardise_out$X_center
      y = standardise_out$y
      y_mean = standardise_out$y_mean
      scale_pen = standardise_out$scale_pen
      rm(standardise_out)
      if (intercept & type=="logistic"){ 
        X_inter = Matrix::cbind2(1,X)
        groups_inter = c(1,groups+1)
        p_inter = num_vars+1
      } 
  }
  
  # -------------------------------------------------------------
  # weights 
  # ------------------------------------------------------------- 
  if (model == "sgs"){
    if (is.null(v_weights) & is.null(w_weights)){
      pens_out = gen_pens(gFDR, vFDR, pen_method, groups, alpha)
      pen_slope_org = pens_out$pen_slope_org
      pen_gslope_org = pens_out$pen_gslope_org
    } else {
      pen_slope_org = v_weights
      pen_gslope_org = w_weights
    }
  } else if (model == "gslope") {
    pen_slope_org = rep(0, num_vars)
    if (is.null(w_weights)){
      pens_out = gen_pens(gFDR, 0.1, pen_method, groups, 0)
      pen_gslope_org = pens_out$pen_gslope_org
    } else {
      pen_gslope_org = w_weights
    } 
    if (!is.decreasing(pen_gslope_org)){
    stop("penalty sequences must be decreasing")
    }
  }

  # weight checks
  if (any(pen_slope_org < 0) | any(pen_gslope_org < 0)){
    stop("penalty sequences must be non-negative")
  }

  # -------------------------------------------------------------
  # path 
  # -------------------------------------------------------------
  if (any(lambda == "path")){
    lambda_path = path_fcn(grad_vec_zero, groups, groupIDs, alpha, pen_gslope_org, pen_slope_org, path_length, min_frac, sqrt(table(groups)))
  } else if (fit_type == "single"){ 
    if (screen){
      lambda_path = path_fcn(grad_vec_zero, groups, groupIDs, alpha, pen_gslope_org, pen_slope_org, 3, 0.2, sqrt(table(groups)))
      if (lambda_path[1] <= lambda){
          warning("input lambda is larger than lambda max. No screening possible")
      }
      lambda_path = c(lambda_path[1],lambda)
    } else {lambda_path = lambda}
      path_length = length(lambda_path)
  } else if (!any(lambda == "path") & length(lambda) > 1){
    lambda_path = lambda
    path_length = length(lambda_path)
  }

  # -------------------------------------------------------------
  # fitting
  # ------------------------------------------------------------- 
  fitting_options = list(num_obs=num_obs, max_iter=max_iter, backtracking=backtracking, max_iter_backtracking = max_iter_backtracking,
                          f=f,f_grad=f_grad,mult_fcn=mult_fcn,crossprod_mat=crossprod_mat,tol=tol,verbose=verbose)
 
  if (screen){ # screening
    screen_out = do.call(screen_strong, c(list(X, y, groups, groupIDs, type, lambda_path*scale_pen, alpha, pen_slope_org, pen_gslope_org, 
                                        X_scale, num_vars, wt, path_length, model, var_screen_fcn, grp_screen_fcn, kkt_fcn), fitting_options))
    if (fit_type == "single"){ # basic rule
      if (model == "gslope"){
        screen_out$screen_set = screen_out$screen_set_grp[[-1]]
        screen_out$kkt_violations = screen_out$kkt_violations_grp[[-1]]
        screen_out$epsilon_set = screen_out$epsilon_set_grp[[-1]]
      } else {
        screen_out$screen_set_var = screen_out$screen_set_var[[-1]]
        screen_out$screen_set_grp = screen_out$screen_set_grp[[-1]]
        screen_out$kkt_violations_var = screen_out$kkt_violations_var[[-1]]
        screen_out$kkt_violations_grp = screen_out$kkt_violations_grp[[-1]]
        screen_out$epsilon_set_var = screen_out$epsilon_set_var[[-1]]
        screen_out$epsilon_set_grp = screen_out$epsilon_set_grp[[-1]]
      }
      screen_out$beta = screen_out$beta[,-1]
      screen_out$group_effects = screen_out$group_effects[,-1]
      screen_out$active_set_var = screen_out$active_set_var[[-1]]
      screen_out$active_set_grp = screen_out$active_set_grp[[-1]]
      screen_out$num_it = screen_out$num_it[-1]
      screen_out$success = screen_out$success[-1]
      screen_out$certificate = screen_out$certificate[-1]
      screen_out$x = screen_out$x[,-1]
      screen_out$z = screen_out$z[,-1]
      screen_out$u = screen_out$u[,-1]
      lambda_path = lambda_path[2]
    } 
    # prepare output
    out = c()
    out$beta = screen_out$beta
    out$group_effects = screen_out$group_effects
    out$selected_var = screen_out$active_set_var
    out$selected_grp = screen_out$active_set_grp
    if (intercept){
      out$beta = apply(out$beta,2,function(x) c(y_mean - sum(X_center*x),x))
    } 
    out$num_it = screen_out$num_it
    out$success = screen_out$success
    out$certificate = screen_out$certificate
    out$x = screen_out$x
    out$z = screen_out$z
    out$u = screen_out$u
    if (model == "gslope"){
      out$screen_set = screen_out$screen_set_grp
      out$epsilon_set = screen_out$epsilon_set_grp
      out$kkt_violations = screen_out$kkt_violations_grp
    } else {
      out$screen_set_var = screen_out$screen_set_var
      out$screen_set_grp = screen_out$screen_set_grp
      out$epsilon_set_var = screen_out$epsilon_set_var
      out$epsilon_set_grp = screen_out$epsilon_set_grp
      out$kkt_violations_var = screen_out$kkt_violations_var
      out$kkt_violations_grp = screen_out$kkt_violations_grp
    }
  } else { # no screening 
  if (fit_type == "single"){
    out = do.call(fit_one, c(list(X,y,groups, groupIDs, type, scale_pen*lambda_path, alpha, intercept, pen_slope_org, pen_gslope_org, rep(0,ncol(X)), rep(0,ncol(X)), X_scale, X_center, y_mean, wt, sqrt(table(groups)) ,model), fitting_options))
    } else {
    out = do.call(fit_path, c(list(X,y,groups, groupIDs, type, scale_pen*lambda_path, alpha, intercept, pen_slope_org, pen_gslope_org, X_scale, X_center, y_mean, wt, sqrt(table(groups)), model, num_vars, num_groups, path_length), fitting_options))
   }
    if (model == "gslope"){
      out$screen_set = "no screening performed"
      out$epsilon_set = "no screening performed"
      out$kkt_violations = "no screening performed"
    } else {
      out$screen_set_var = "no screening performed"
      out$screen_set_grp = "no screening performed"
      out$epsilon_set_var = "no screening performed"
      out$epsilon_set_grp = "no screening performed"
      out$kkt_violations_var = "no screening performed"
      out$kkt_violations_grp = "no screening performed"
    }
  }

  # -------------------------------------------------------------
  # prepare output
  # ------------------------------------------------------------- 
  if (model == "gslope"){
    out$pen_slope = c()
    out$pen_gslope = pen_gslope_org
  } else {
    out$pen_slope = pen_slope_org
    out$pen_gslope = pen_gslope_org
  }
  out$screen = screen
  out$type = type
  out$standardise = standardise
  out$intercept = intercept 
  out$lambda = lambda_path
  out$beta = as(out$beta,"CsparseMatrix")
  out$group_effects = as(out$group_effects,"CsparseMatrix")
  rownames(out$group_effects) = paste0("G", 1:num_groups)
  class(out) <- model
  return(out)
}

general_fit_cv = function(X, y, groups, model, path_fcn, type, lambda, path_length, nfolds, alpha, vFDR, gFDR, pen_method, backtracking, max_iter, max_iter_backtracking, tol, min_frac, standardise, intercept, v_weights, w_weights, error_criteria, screen, verbose, gamma_1, gamma_2){
  num_vars = ncol(X)
  num_obs = nrow(X)
  num_groups = length(unique(groups))
  groupIDs = getGroupID(groups) 
  
  # type of model
  if (type == "linear"){ 
    f = mse_loss
    f_grad = mse_grad
  } else if (type == "logistic"){
    f = log_loss
    f_grad = log_grad
  } else {stop("loss function not supported")} 

  # -------------------------------------------------------------
  # weights
  # ------------------------------------------------------------- 
  if (model == "sgs"){
    if (is.null(v_weights) & is.null(w_weights)){
      pens_out = gen_pens(gFDR, vFDR, pen_method, groups, alpha)
      pen_slope_org = pens_out$pen_slope_org
      pen_gslope_org = pens_out$pen_gslope_org
    } else {
      pen_slope_org = v_weights
      pen_gslope_org = w_weights
    }
    if (!is.decreasing(pen_slope_org) | !is.decreasing(pen_gslope_org)){
      stop("penalty sequences must be decreasing")
    }
  } else if (model == "gslope") {
    if (is.null(w_weights)){
      pens_out = gen_pens(gFDR, 0.1, pen_method, groups, 0)
      pen_slope_org = rep(0, num_groups)
      pen_gslope_org = pens_out$pen_gslope_org
    } else {
      pen_slope_org = rep(0, num_groups)
      pen_gslope_org = w_weights
    } 
    if (!is.decreasing(pen_gslope_org)){
    stop("penalty sequences must be decreasing")
    }
  }

  if (any(pen_slope_org < 0) | any(pen_gslope_org < 0)){
    stop("penalty sequences must be non-negative")
  }

  # -------------------------------------------------------------
  # path
  # ------------------------------------------------------------- 
  if (any(lambda == "path")){
    if (intercept){
       grad_vec_zero = arma_mv(t(Matrix::cbind2(1,X)), f_grad(y=y, arma_mv(Matrix::cbind2(1,X), rep(0,num_vars+1)), num_obs = num_obs))[-1]
    } else {grad_vec_zero = arma_mv(Matrix::t(X),f_grad(y=y, arma_mv(X, rep(0,num_vars)), num_obs = num_obs))}
    lambda_path = path_fcn(grad_vec_zero, groups, groupIDs, alpha, pen_gslope_org, pen_slope_org, path_length, min_frac, sqrt(table(groups)))
  } else {
    path_length = length(lambda)
    lambda_path = lambda
  }

  # -------------------------------------------------------------
  # fitting final models
  # ------------------------------------------------------------- 
  # initialise CV variable for storing results
  all_data = data.frame(y,X)
  folds = createFolds(y, k = nfolds, list=TRUE)
  all_errors = matrix(0,nrow=path_length,ncol=nfolds)
  output_errors = data.frame(lambda=lambda_path,error_criteria=rep(0,path_length), num_non_zero = rep(0,path_length))

  if (model == "sgs"){
    lambda_model = fit_sgs(X=X, y=y, groups=groups, type=type, lambda=lambda_path, path_length=path_length, alpha=alpha, max_iter = max_iter, backtracking = backtracking, max_iter_backtracking = max_iter_backtracking, tol = tol, standardise=standardise, intercept=intercept, screen = screen, verbose = FALSE, v_weights = pen_slope_org, w_weights = pen_gslope_org)
  } else if (model == "gslope"){
    lambda_model = fit_gslope(X=X, y=y, groups=groups, type=type, lambda=lambda_path, path_length=path_length,max_iter=max_iter, backtracking=backtracking, max_iter_backtracking=max_iter_backtracking, tol=tol, standardise=standardise, intercept=intercept, w_weights=pen_gslope_org, screen=screen, verbose=FALSE)
  }

  # -------------------------------------------------------------
  # cv loop
  # ------------------------------------------------------------- 
  for (fold_id in 1:nfolds){
    # Create training design matrix and target data, leaving one out each time
    Train = all_data[as.numeric(unlist(folds[-fold_id])),]
    Train_y = Train$y
    Train_X = Train[,-1]

    # Create testing design matrix and target data
    Test = all_data[as.numeric(unlist(folds[fold_id])),]
    Test_y = Test$y
    Test_X = X[as.numeric(unlist(folds[fold_id])),]
  
    # Fit Model
    if (model == "sgs"){
      cv_model = fit_sgs(X=as.matrix(Train_X), y=as.matrix(Train_y), groups=groups, type=type, lambda=lambda_path, alpha=alpha, max_iter=max_iter, backtracking=backtracking, max_iter_backtracking=max_iter_backtracking, tol=tol, standardise=standardise, intercept=intercept,
         v_weights= pen_slope_org, w_weights = pen_gslope_org, screen = screen, verbose = FALSE)
    } else if (model == "gslope"){
      cv_model = fit_gslope(X=as.matrix(Train_X), y=as.matrix(Train_y), groups=groups, type=type, lambda=lambda_path, max_iter=max_iter, backtracking=backtracking, max_iter_backtracking=max_iter_backtracking, tol=tol, standardise=standardise, intercept=intercept,
         w_weights = pen_gslope_org, screen = screen, verbose = FALSE)
    }
     
    # Error
    if (type=="linear"){
      if (intercept){
        if (error_criteria == "mse"){
          error_val = apply(cv_model$beta,2,function(x) sum((Test_y-arma_mv(Matrix::cbind2(1,Test_X),as.vector(x)))^2))}
          else if (error_criteria == "mae") {error_val = apply(cv_model$beta,2,function(x) sum(abs(Test_y-arma_mv(Matrix::cbind2(1,Test_X),as.vector(x)))))} else{stop("not a valid criteria")}
      } else {
        if (error_criteria == "mse"){
          error_val = apply(cv_model$beta, 2, function(x) sum((Test_y-arma_mv(Test_X,as.vector(x)))^2))
          } else if (error_criteria =="mae"){error_val = apply(cv_model$beta, 2, function(x) sum(abs(Test_y-arma_mv(Test_X,as.vector(x)))))} else {stop("not a valid criteria")}
      }
      } else if (type=="logistic"){
        if (intercept){
          error_val = apply(cv_model$beta,2,function(x) 1-sum(ifelse(sigmoid(arma_mv(Matrix::cbind2(1,Test_X),as.vector(x)))>=0.5,1,0) == Test_y)/length(Test_y))
        } else {
          error_val = apply(cv_model$beta,2,function(x) 1-sum(ifelse(sigmoid(arma_mv(Test_X,as.vector(x)))>=0.5,1,0) == Test_y)/length(Test_y))
        }
      }

    all_errors[, fold_id] = error_val
    if (verbose == TRUE){print(paste0("Fold ", fold_id,"/",nfolds, " done. Error: ", mean(error_val)))}
  }

  # -------------------------------------------------------------
  # prepare output
  # ------------------------------------------------------------- 
  output_errors$error_criteria = apply(all_errors,1,mean)
  output_errors$num_non_zero = unlist(lapply(lambda_model$selected_var, length))

  # Pick best lambda - 1se
  error_se = apply(all_errors,1,sd)/sqrt(nfolds)
  error_se_lambda_min = error_se[which.min(output_errors$error_criteria)]
  best_lambda = max(lambda_path[output_errors$error_criteria < min(output_errors$error_criteria) + error_se_lambda_min])

  # prepare output for best model
  output_model = c()
  best_lambda_id = match(best_lambda, lambda_path)
  output_model$beta = lambda_model$beta[,best_lambda_id]
  output_model$group_effects = lambda_model$group_effects[,best_lambda_id]
  output_model$selected_var = lambda_model$selected_var[[best_lambda_id]]
  output_model$selected_grp = lambda_model$selected_grp[[best_lambda_id]]
  output_model$pen_slope = lambda_model$pen_slope
  output_model$pen_gslope = lambda_model$pen_gslope
  output_model$lambda = lambda_model$lambda[best_lambda_id]
  output_model$type = lambda_model$type
  output_model$intercept = lambda_model$intercept
  output_model$standardise = lambda_model$standardise
  output_model$num_it = lambda_model$num_it[best_lambda_id]
  output_model$success = lambda_model$success[best_lambda_id]
  output_model$certificate = lambda_model$certificate[best_lambda_id]
  output_model$x = lambda_model$x[,best_lambda_id]
  output_model$z = lambda_model$z[,best_lambda_id]
  output_model$u = lambda_model$u[,best_lambda_id]
  output_model$screen = lambda_model$screen
  if (model == "gslope"){
    if (lambda_model$screen){
      output_model$screen_set = lambda_model$screen_set[[best_lambda_id]]
      output_model$epsilon_set = lambda_model$epsilon_set[[best_lambda_id]]
      output_model$kkt_violations = lambda_model$kkt_violations[[best_lambda_id]]
    } else {
      output_model$screen_set = lambda_model$screen_set[[1]]
      output_model$epsilon_set = lambda_model$epsilon_set[[1]]
      output_model$kkt_violations = lambda_model$kkt_violations[[1]] 
    }
  } else {
    if (lambda_model$screen){
      output_model$screen_set_var = lambda_model$screen_set_var[[best_lambda_id]]
      output_model$screen_set_grp = lambda_model$screen_set_grp[[best_lambda_id]]
      output_model$epsilon_set_var = lambda_model$epsilon_set_var[[best_lambda_id]]
      output_model$epsilon_set_grp = lambda_model$epsilon_set_grp[[best_lambda_id]]
      output_model$kkt_violations_var = lambda_model$kkt_violations_var[[best_lambda_id]]
      output_model$kkt_violations_grp = lambda_model$kkt_violations_grp[[best_lambda_id]]
    } else {
      output_model$screen_set_var = lambda_model$screen_set_var[[1]]
      output_model$screen_set_grp = lambda_model$screen_set_grp[[1]]
      output_model$epsilon_set_var = lambda_model$epsilon_set_var[[1]]
      output_model$epsilon_set_grp = lambda_model$epsilon_set_grp[[1]]
      output_model$kkt_violations_var = lambda_model$kkt_violations_var[[1]]
      output_model$kkt_violations_grp = lambda_model$kkt_violations_grp[[1]] 
    }
  }
  out = c()
  out$all_models = lambda_model
  out$fit = output_model
  out$best_lambda = best_lambda
  out$best_lambda_id = best_lambda_id
  out$errors = output_errors
  class(out) <- paste0(model,"_cv")
  return(out)
}