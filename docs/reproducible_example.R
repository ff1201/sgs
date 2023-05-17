## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,message=FALSE ,warning=FALSE)

## -----------------------------------------------------------------------------
library(sgs)
groups = c(rep(1:20, each=3),
           rep(21:40, each=4),
           rep(41:60, each=5),
           rep(61:80, each=6),
           rep(81:100, each=7))

data = generate_toy_data(p=500, n=400, groups = groups, seed_id=3)

## -----------------------------------------------------------------------------
model = fit_sgs(X = data$X, y = data$y, groups = groups, type="linear", lambda = 1, alpha=0.95, vFDR=0.1, gFDR=0.1, standardise = "l2", intercept = TRUE, verbose=FALSE)

## -----------------------------------------------------------------------------
model$beta[model$selected_var+1,]
model$group.effects[model$selected_group,]
model$selected_var
model$selected_group

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
fdr_sensitivity(fitted_ids = model$selected_var, true_ids = data$true_var_id, num_coef = 500)
fdr_sensitivity(fitted_ids = model$selected_group, true_ids = data$true_grp_id, num_coef = 100)

## -----------------------------------------------------------------------------
cv_model = fit_sgs_cv(X = data$X, y = data$y, groups=groups, type = "linear", nlambda = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, standardise="l2",intercept=TRUE,verbose=TRUE)

## -----------------------------------------------------------------------------
print(cv_model)

## -----------------------------------------------------------------------------
fdr_sensitivity(fitted_ids = cv_model$fit$selected_var, true_ids = data$true_var_id, num_coef = 500)
fdr_sensitivity(fitted_ids = cv_model$fit$selected_group, true_ids = data$true_grp_id, num_coef = 100)

## -----------------------------------------------------------------------------
plot(cv_model,how_many = 10)

## -----------------------------------------------------------------------------
predict(model,data$X,type="linear")[1:5]

## -----------------------------------------------------------------------------
sigmoid = function(x) {
   1 / (1 + exp(-x))
}
y = ifelse(sigmoid(data$X %*% data$true_beta + rnorm(400))>0.5,1,0)
train_y = y[1:350] 
test_y = y[351:400]
train_X = data$X[1:350,] 
test_X = data$X[351:400,]

## -----------------------------------------------------------------------------
cv_model = fit_sgs_cv(X = train_X, y = train_y, groups=groups, type = "logistic", nlambda = 20, nfolds=10, alpha = 0.95, vFDR = 0.1, gFDR = 0.1, min_frac = 0.05, standardise="l2",intercept=FALSE,verbose=TRUE)

## -----------------------------------------------------------------------------
predictions = predict(cv_model$fit,test_X,type="logistic")

## -----------------------------------------------------------------------------
predictions$response[1:5]
predictions$class[1:5]
sum(predictions$class == test_y)/length(test_y)

