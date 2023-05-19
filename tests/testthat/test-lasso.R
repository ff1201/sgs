test_that("solution reduces to lasso when alpha=1 and constant weights, with no intercept or standardisation", {
  library(glmnet)
  n = 50
  p = 100
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso <- glmnet(X, y, lambda = lambda, standardize = FALSE,family="gaussian",intercept=FALSE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=1, vFDR=0.1, gFDR=0.1,standardise="none",intercept=FALSE,w_weights = rep(0,p),v_weights = rep(lambda,p),tol=1e-5)
    
  lasso_cost = sgs_convex_opt(X=X,y=y,beta= as.matrix(lasso$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(as.matrix(lasso$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(lasso_cost,sgs_cost, tol=1e-3)
})

test_that("solution reduces to lasso when alpha=1 and constant weights, with intercept", {
  library(glmnet)
  n = 50
  p = 100
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso <- glmnet(X, y, lambda = lambda, standardize = FALSE,family="gaussian",intercept=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=1, vFDR=0.1, gFDR=0.1,standardise="none",intercept=TRUE,w_weights = rep(0,p),v_weights = rep(lambda,p),tol=1e-5)
    
  lasso_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=TRUE)

  expect_equivalent(as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(lasso_cost,sgs_cost, tol=1e-3)
})


test_that("solution reduces to lasso when alpha=1 and constant weights, using standardisation but no intercept", {
  library(glmnet)
  n = 50
  p = 100
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  X = scale(X,center=TRUE,scale=FALSE) # intercept=TRUE centers X in glmnet
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso <- glmnet(X, y, lambda = lambda, standardize = TRUE,family="gaussian",intercept=FALSE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=1, vFDR=0.1, gFDR=0.1,standardise="sd",intercept=FALSE,w_weights = rep(0,p),v_weights = rep(lambda,p),tol=1e-5)
    
  lasso_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(lasso$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=FALSE)

  expect_equivalent(as.matrix(lasso$beta),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(lasso_cost,sgs_cost, tol=1e-3)
})

test_that("solution reduces to lasso when alpha=1 and constant weights, using standardisation and intercept", { # sd off by a very small amount
  library(glmnet)
  n = 50
  p = 100
  data= generate_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda = 0.8
  groups = 1:p
  lasso <- glmnet(X, y, lambda = lambda, standardize = TRUE,family="gaussian",intercept=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=1, alpha=1, vFDR=0.1, gFDR=0.1,standardise="sd",intercept=TRUE,w_weights = rep(0,p),v_weights = rep(lambda,p),tol=1e-5)
    
  lasso_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$pen_gslope,slope_seq=sgs$pen_slope,groups=groups,intercept=TRUE)

  expect_equivalent(as.matrix(c(as.matrix(lasso$a0), as.matrix(lasso$beta))),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(lasso_cost,sgs_cost, tol=1e-3)
})