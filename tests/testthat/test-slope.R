test_that("solution reduces to slope when alpha=1, with no intercept or standardisation", {
  n = 50
  p = 100
  data= gen_toy_data(p=p,n=n,rho = 0,seed_id = 3,grouped = FALSE,var_sparsity=0.2,orthogonal = FALSE)
  X <- data$X
  y <- data$y
  lambda=0.8
  groups=1:p

  coef_ref <- c(
    0, 0, 6.35351249775198, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.924615812626428,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.593341879455522, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.0574674799067196, 0, 0, 0, 0, 0, 0, 0, 0
  )

  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="none")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    -0.591904300904075, 0, 0, 6.15570369828865, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -0.855516332321258, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.597393440042509,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0.147583863629151, 0, 0, 0, 0, 0, 0,
    0, 0
  )

  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="none")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    0, 0, 6.62193929302955, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.794831421813014,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.446180314946519, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.000446641410816172, 0, 0, 0, 0, 0, 0, 0, 0
  )

  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="sd")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    0, 0, 7.08023698435155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.32172932399327, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.72816894344401, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0.281886263999637, 0, 0, 0, 0.421166293233546, 0, 0, 0, 0, 0, 
    0, 0, 0
  )

  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="l1")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(coef_ref,
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

  slope = SLOPE::SLOPE(X, y, family = "gaussian", alpha = lambda, q=0.1,intercept=FALSE,scale="l2",center=TRUE)

  coef_ref <- c(
    0, 0, 6.62199573208823, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.794796341616357,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.446206778892696, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.000340786752155744, 0, 0, 0, 0, 0, 0, 0, 0
  )
 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=FALSE,standardise="l2")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=FALSE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    -0.443697677962616, 0, 0, 6.621939698726, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.794830608447868,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.446179964202237, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.000444025886934943, 0, 0, 0, 0, 0, 0, 0, 0
  )
 
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="sd")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    -0.237185741459661, 0, 0, 7.08023698435156, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -1.32172932399327, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.728168943444005,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0.281886263999636, 0, 0, 0, 0.421166293233551, 0,
    0, 0, 0, 0, 0, 0, 0
  )
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="l1")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(coef_ref,
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

  coef_ref <- c(
    -0.443672029601005, 0, 0, 6.62199573208823, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    -0.794796341616357, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.446206778892693,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0.00034078675215286, 0, 0, 0, 0, 0, 0,
    0, 0
  )
  
  sgs = fit_sgs(X=X,y=y, groups=groups, type="linear", lambda=lambda, alpha=1, vFDR=0.1, gFDR=0.1,intercept=TRUE,standardise="l2")
  
  slope_cost = sgs_convex_opt(X=X,y=y,beta=coef_ref,num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)
  sgs_cost = sgs_convex_opt(X=X,y=y,beta=as.matrix(sgs$beta),num_obs=n,gslope_seq=sgs$w_weights,slope_seq=sgs$v_weights,groups=groups, intercept=TRUE)

  expect_equivalent(coef_ref,
    as.matrix(sgs$beta),
    tol = 1e-3
  )

  expect_equivalent(slope_cost, sgs_cost, tol=1e-3)
})
