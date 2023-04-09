source("fit_sgs.R")

fit_sgs_cv = function(X,y,groups,pen_method=1, type = "linear", nlambda = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, backtracking = 0.7, max_iter = 5000, max_iter_backtracking = 100, tol = 1e-5, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=FALSE, v_weights=NULL,w_weights=NULL,error_criteria = "mse",max_lambda = NULL){
  # generate lambdas
  p = ncol(X)
  num_obs = nrow(X)
  num_groups = length(unique(groups))
  group_ids = getGroupID(groups) 
  len_each_grp = sapply(group_ids, length)
  wt_per_grp = sqrt(len_each_grp)
  
  # weights
  if (is.null(v_weights) & is.null(w_weights)){
    pens_out = generate_penalties(gFDR, vFDR, pen_method,num_groups,wt_per_grp, len_each_grp, p,alpha)
    pen_slope_org = pens_out$pen_slope_org
    pen_gslope_org = pens_out$pen_gslope_org
  } else {
    pen_slope_org = v_weights
    pen_gslope_org = w_weights
  }

  # standardise data
  standardise_out = standardise_sgs(X=X,y=y,standardise=standardise,intercept=intercept,num_obs=num_obs)
  X_path = standardise_out$X
  y_path = standardise_out$y
  
  # calculate path
  if (is.null(max_lambda)){
    lambdas = (1/standardise_out$scale_pen)*generate_lambda_path(X=as.matrix(X_path),y=as.matrix(y_path),groups=groups,alpha=alpha,min_frac=min_frac,nlambda=nlambda, v_weights=pen_slope_org, w_weights=pen_gslope_org,group.sizes = len_each_grp)
  } else {
    min_lambda = min_frac*max_lambda
    lambdas = exp(seq(log(max_lambda),log(min_lambda), (log(min_lambda) - log(max_lambda))/(nlambda-1))) 
  }

  # initialise CV variable for storing results
  all_data = data.frame(y,X)
  set.seed(5)
  folds = createFolds(y, k = nfolds, list=TRUE)
  all_errors = matrix(0,nrow=nlambda,ncol=nfolds)
  output_errors = data.frame(lambda=lambdas,error_criteria=rep(0,nlambda), num_non_zero = rep(0,nlambda))
  ## warm starts
  initial_x = rep(0,p)
  initial_u = rep(0,p)

  list_of_models = list()
  
  for (lambda_id in 1:nlambda){
    num_its = 0
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
      model = fit_sgs(X=as.matrix(Train_X), y=as.matrix(Train_y),groups=groups, type=type, lambda=lambdas[lambda_id], alpha=alpha, max_iter = max_iter, backtracking = backtracking, max_iter_backtracking = max_iter_backtracking, tol = tol, standardise=standardise, intercept=intercept,
        x0 = initial_x, u = initial_u, v_weights= pen_slope_org, w_weights = pen_gslope_org)

      # Error
      if (type=="linear"){
      if (intercept){
        if (error_criteria == "mse"){
          error_val = sum((Test_y-arma_mm(cbind(1,Test_X),as.vector(model$beta)))^2)}
          else if (error_criteria == "mae") {error_val = sum(abs(Test_y-arma_mm(cbind(1,Test_X),as.vector(model$beta))))} else{stop("not a valid criteria")}
      } else {
        if (error_criteria == "mse"){
          error_val = sum((Test_y-arma_mm(Test_X,as.vector(model$beta)))^2)
          } else if (error_criteria =="mae"){error_val = sum(abs(Test_y-arma_mm(Test_X,as.vector(model$beta))))} else {stop("not a valid criteria")}
      }
      } else if (type=="logistic"){
        if (intercept){
          error_val = 1-sum(ifelse(sigmoid(arma_mm(cbind(1,Test_X),as.vector(model$beta)))>=0.5,1,0) == Test_y)/length(Test_y)
        } else {
          error_val = 1-sum(ifelse(sigmoid(arma_mm(Test_X,as.vector(model$beta)))>=0.5,1,0) == Test_y)/length(Test_y)
        }
      }
      all_errors[lambda_id, fold_id] = error_val
      num_its = num_its + model$num_it
  }
      lambda_model = fit_sgs(X=X, y=y,groups=groups, type=type, lambda=lambdas[lambda_id], alpha=alpha, max_iter = max_iter, backtracking = backtracking, max_iter_backtracking = max_iter_backtracking, tol = tol, standardise=standardise, intercept=intercept,
      x0 = initial_x, u = initial_u, v_weights= pen_slope_org, w_weights = pen_gslope_org)
      list_of_models[[lambda_id]] = lambda_model
      # warm starts
      initial_x = lambda_model$x
      initial_u = lambda_model$u
      output_errors$error_criteria[lambda_id] = mean(all_errors[lambda_id,])
      output_errors$num_non_zero[lambda_id] = length(lambda_model$selected_var)
      if (verbose == TRUE){
        if (type == "linear"){print(paste0("Lambda ", lambda_id,"/",nlambda, " done. Lambda: ", round(lambdas[lambda_id],4), ". Number of non-zero: ",length(lambda_model$selected_var),". Error: ", output_errors$error_criteria[lambda_id], ". Avg iter: ", floor(num_its/nfolds)))}
        else if (type == "logistic"){
          print(paste0("Lambda ", lambda_id,"/",nlambda, " done. Lambda: ", round(lambdas[lambda_id],4), ". Number of non-zero: ",length(lambda_model$selected_var),". Misclassification error: ", output_errors$error_criteria[lambda_id], ". Avg iter: ", floor(num_its/nfolds)))
          } }
  }

  # Pick best lambda - 1se
  error_se = apply(all_errors,1,sd)/sqrt(nfolds)
  error_se_lambda_min = error_se[which.min(output_errors$error_criteria)]
  best_lambda = max(lambdas[output_errors$error_criteria < min(output_errors$error_criteria) + error_se_lambda_min])
 
  model = list_of_models[[match(best_lambda, lambdas)]]
  out = c()
  out$all_models =list_of_models
  out$fit = model
  out$best_lambda = best_lambda
  out$best_lambda_id = match(best_lambda, lambdas)
  out$errors = output_errors
  out$type = type
  class(out) <- "SGS_CV"
  return(out)
}

print.SGS_CV = function(out, digits = max(3, getOption("digits") - 3), ...){ # fix this
  cat("\n regression type: ", out$type, "\n\n")
  print(cbind(lambda = out$errors$lambda, error = out$errors$error_criteria, estimated_non_zero = out$errors$num_non_zero))
}