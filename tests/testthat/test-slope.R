test_that("solution reduces to slope when alpha=1, with no intercept or standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=FALSE,solver="admm",screen=FALSE,scale="none",center=FALSE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with intercept but no standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=TRUE,solver="admm",screen=FALSE,scale="none",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with no intercept but sd standardisation", { # again, sd off by a small amount
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=FALSE,solver="admm",screen=FALSE,scale="sd",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="sd")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with no intercept but l1 standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=FALSE,solver="admm",screen=FALSE,scale="l1",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="l1")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with no intercept but l2 standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=FALSE,solver="admm",screen=FALSE,scale="l2",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="l2")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with intercept and sd standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=TRUE,solver="admm",screen=FALSE,scale="sd",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="sd")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with intercept and l1 standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=TRUE,solver="admm",screen=FALSE,scale="l1",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="l1")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})

test_that("solution reduces to slope when alpha=1, with intercept and l2 standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=TRUE,solver="admm",screen=FALSE,scale="l2",center=TRUE)
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="l2")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=slope$coefficients,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(as.matrix(slope$coefficients),
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})